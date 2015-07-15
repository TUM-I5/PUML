/**
 * @file
 *  This file is part of PUML
 *
 *  For conditions of distribution and use, please see the copyright
 *  notice in the file 'COPYING' at the root directory of this package
 *  and the copyright notice at https://github.com/TUM-I5/PUML
 *
 * @copyright 2013 Technische Universitaet Muenchen
 * @author Sebastian Rettenberger <rettenbs@in.tum.de>
 */

#include <mpi.h>

#include <algorithm>
#include <cstring>
#include <map>
#include <string>
#include <vector>

#include <parmetis.h>

#include <netcdf_par.h>
#include <netcdf.h>

#include <apfMesh2.h>
#include <apfNumbering.h>
#include <apfZoltan.h>
#include <maMesh.h>

#include "utils/args.h"
#include "utils/logger.h"

#include "input/SerialMeshFile.h"
#include "input/ApfNative.h"
#ifdef USE_SIMMOD
#include "input/SimModSuite.h"
#endif // USE_SIMMOD
#include "meshreader/ParallelGambitReader.h"
#include "meshreader/ParallelFidapReader.h"

const static unsigned int FACE2INTERNAL[] = {0, 1, 3, 2};

struct MPINeighborElement {
	/** Local element */
	apf::MeshEntity* element;
	/**	Global id of the local element */
	long localId;
	/** Side of the local element */
	int localSide;
	/** Global number neighbor element */
	long neighborId;
	/** Side of the neighbor element */
	int neighborSide;
};

static void checkNcError(int error)
{
	if (error != NC_NOERR)
		logError() << "Error while writing netCDF file:" << nc_strerror(error);
}

bool compareLocalMPINeighbor(const MPINeighborElement &elem1, const MPINeighborElement &elem2)
{
	return (elem1.localId < elem2.localId)
			|| (elem1.localId == elem2.localId && FACE2INTERNAL[elem1.localSide] < FACE2INTERNAL[elem2.localSide]);
}

bool compareRemoteMPINeighbor(const MPINeighborElement &elem1, const MPINeighborElement &elem2)
{
	return (elem1.neighborId < elem2.neighborId)
			|| (elem1.neighborId == elem2.neighborId && FACE2INTERNAL[elem1.neighborSide] < FACE2INTERNAL[elem2.neighborSide]);
}

/** Maps from element + face to orientation */
const static int FACE2ORIENTATION[4][4] = {
		{0, 2, 1, -1},
		{0, 1, -1, 2},
		{-1, 0, 1, 2},
		{0, -1, 2, 1}
};

static int getOrientation(const long* localVertNum, int localFace, long firstNbVert)
{
	unsigned int off = std::find(localVertNum, localVertNum+4, firstNbVert)-localVertNum;
	assert(off >= 0 && off <= 3);
	int o = FACE2ORIENTATION[localFace][off];
	assert(o >= 0);
	return o;
}

static apf::MeshEntity* getFaceElemOppositeElem(apf::Mesh* m,
		apf::MeshEntity* face, apf::MeshEntity* elem)
{
	apf::Up up;
	m->getUp(face, up);
	if (up.n <= 1)
		return 0L;
	if (up.e[0] == elem)
		return up.e[1];
	return up.e[0];
}

/**
 * @return The number of the first vertex for face <code>face</code>
 *  of the element <code>element</code>
 */
static long getFirstNumOfFace(apf::GlobalNumbering* num,
		apf::MeshEntity* element, int face)
{
	apf::NewArray<long> vn;
	apf::getElementNumbers(num, element, vn);

	return vn[face == 2 ? 1 : 0];
}

const static int METIS_RANDOM_SEED = 42;

int main(int argc, char* argv[])
{
	int rank = 0;
	int processes = 1;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &processes);

	PCU_Comm_Init();

	// Parse command line arguments
	utils::Args args;
	const char* source[] = {"gambit", "fidap", "apf", "simmodsuite"};
	args.addEnumOption("source", source, 's', "Mesh source (default: gambit)", false);
	args.addOption("dump", 'd', "Dump APF mesh before partitioning it",
			utils::Args::Required, false);
	args.addOption("model", 0, "Dump/Load a specific model file",
			utils::Args::Required, false);
	args.addOption("vtk", 0, "Dump mesh to VTK files",
			utils::Args::Required, false);
	args.addOption("license", 'l', "License file (only used by SimModSuite)",
			utils::Args::Required, false);
	args.addOption("cad", 'c', "CAD file (only used by SimModSuite)",
			utils::Args::Required, false);
	args.addOption("mesh", 0, "Mesh attributes name (only used by SimModSuite, default: \"mesh\")",
			utils::Args::Required, false);
	args.addOption("analysis", 0, "Analysis attributes name (only used by SimModSuite, default: \"analysis\")",
			utils::Args::Required, false);
	args.addOption("stl", 0, "Use STL-input with hard-coded parameters",
			utils::Args::No, false);
	const char* forces[] = {"0", "1", "2"};
	args.addEnumOption("enforce-size", forces, 0, "Enforce mesh size (only used by SimModSuite, default: 0)", false);
	args.addAdditionalOption("input", "Input file (mesh or model)");
	args.addAdditionalOption("partition", "Number of partitions");
	args.addAdditionalOption("output", "Output parallel unstructured mesh file", false);

	if (args.parse(argc, argv, rank == 0) != utils::Args::Success)
		return 1;

	const char* inputFile = args.getAdditionalArgument<const char*>("input");

	unsigned int nPartitions = args.getAdditionalArgument<unsigned int>("partition");
	if (nPartitions == 0)
		logError() << "Partitions created must be greater than zero";

	std::string outputFile;
	if (args.isSetAdditional("output")) {
		outputFile = args.getAdditionalArgument<std::string>("output");
	} else {
		// Compute default output filename
		outputFile = inputFile;
		size_t dotPos = outputFile.find_last_of('.');
		if (dotPos != std::string::npos)
			outputFile.erase(dotPos);
		outputFile.append(".nc.pum");
	}

	MeshInput* meshInput = 0L;
	apf::Mesh2* mesh = 0L;
	switch (args.getArgument<int>("source", 0)) {
	case 0:
		logInfo(rank) << "Using Gambit mesh";
		meshInput = new SerialMeshFile<ParallelGambitReader>(inputFile);
		break;
	case 1:
		logInfo(rank) << "Using Fidap mesh";
		meshInput = new SerialMeshFile<ParallelFidapReader>(inputFile);
		break;
	case 2:
		logInfo(rank) << "Using APF native format";
		meshInput = new ApfNative(inputFile,
				args.getArgument<const char*>("model", 0L));
		break;
	case 3:
#ifdef USE_SIMMOD
		logInfo(rank) << "Using SimModSuite";

		meshInput = new SimModSuite(inputFile,
				args.getArgument<const char*>("cad", 0L),
				args.getArgument<const char*>("license", 0L),
				args.getArgument<const char*>("mesh", "mesh"),
				args.getArgument<const char*>("analysis", "analysis"),
				args.getArgument<int>("enforce-size", 0),
				args.isSet("stl"));

#else // USE_SIMMOD
		logError() << "SimModSuite is not supported in this version";
#endif // USE_SIMMOD
		break;
	default:
		logError() << "Unknown source";
	}

	mesh = meshInput->getMesh();

	// Check mesh
	if (alignMdsMatches(mesh))
		logWarning() << "Fixed misaligned matches";
	mesh->verify();

	// Dump mesh for later usage?
	const char* dumpFile = args.getArgument<const char*>("dump", 0L);
	if (dumpFile) {
		logInfo(PCU_Comm_Self()) << "Writing native APF mesh";
		mesh->writeNative(dumpFile);

		const char* modelFile = args.getArgument<const char*>("model", 0L);
		if (modelFile)
			gmi_write_dmg(mesh->getModel(), modelFile);
	}

	// Dump VTK mesh
	const char* vtkPrefix = args.getArgument<const char*>("vtk", 0L);
	if (vtkPrefix) {
		logInfo(PCU_Comm_Self()) << "Writing VTK mesh";
		apf::writeVtkFiles(vtkPrefix, mesh);
	}

	// Get local size
	unsigned int nLocalElements = apf::countOwned(mesh, 3);

	// Get global size
	unsigned int nElements;
	MPI_Allreduce(&nLocalElements, &nElements, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
	logInfo(rank) << "Mesh size:" << nElements;

	// Compute min insphere radius
	double min = std::numeric_limits<double>::max();
	apf::MeshIterator* it = mesh->begin(3);
	while (apf::MeshEntity* element = mesh->iterate(it)) {
		min = std::min(min, ma::getInsphere(mesh, element));
	}
	mesh->end(it);
	MPI_Reduce((rank == 0 ? MPI_IN_PLACE : &min),
			&min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	logInfo(rank) << "Minimum insphere found:" << min;

	// Get element mapping
	unsigned int* elemStart = new unsigned int[processes];
	MPI_Scan(&nLocalElements, &elemStart[rank], 2, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
	elemStart[rank] -= nLocalElements;
	MPI_Allgather(MPI_IN_PLACE, 1, MPI_UNSIGNED, elemStart, 1, MPI_UNSIGNED,
    		MPI_COMM_WORLD);

	// Create dual graph
	logInfo(rank) << "Creating dual graph";
	int* dualGraph = apf::getElementToElement(mesh);

	// Create the partitions
	logInfo(rank) << "Creating partitions with METIS";

	// TODO wrap this code part so we can switch to different libraries
	idx_t* elemDist = new idx_t[processes+1];
	for (int i = 0; i < processes; i++)
		elemDist[i] = elemStart[i];
	elemDist[processes] = nElements;

	idx_t* xadj = new idx_t[nLocalElements+1];
	xadj[0] = 0;
	for (unsigned int i = 0; i < nLocalElements; i++) {
		int neighbors = 0;
		for (unsigned int j = 0; j < 4; j++)
			if (dualGraph[i*4 + j] >= 0)
				neighbors++;
		xadj[i+1] = xadj[i] +  neighbors;
	}

	idx_t* adjncy = new idx_t[xadj[nLocalElements]];
	unsigned int pos = 0;
	for (unsigned int i = 0; i < nLocalElements*4; i++) {
		if (dualGraph[i] >= 0) {
			adjncy[pos] = dualGraph[i];
			pos++;
		}
	}

	delete [] dualGraph;

	idx_t wgtflag = 0;
	idx_t numflag = 0;
	idx_t ncon = 1;
	idx_t nparts = nPartitions;

	real_t* tpwgts = new real_t[nPartitions];
	for (unsigned int i = 0; i < nPartitions; i++)
		tpwgts[i] = 1./nPartitions;

	real_t ubev = 1.01;
	idx_t options[3] = {1, 0, METIS_RANDOM_SEED};

	idx_t edgecut;
	idx_t* part = new idx_t[nLocalElements];

	MPI_Comm commWorld = MPI_COMM_WORLD;
	if (ParMETIS_V3_PartKway(elemDist, xadj, adjncy, 0L, 0L, &wgtflag, &numflag, &ncon,
			&nparts, tpwgts, &ubev, options, &edgecut, part, &commWorld) != METIS_OK)
		logError() << "Could not create partitions";

	delete [] elemDist;
	delete [] xadj;
	delete [] adjncy;
	delete [] tpwgts;

	// Compute the number of partitions for each rank
	unsigned int nMaxLocalPart = (nPartitions + processes - 1) / processes;
	unsigned int nLocalPart = nMaxLocalPart;
	if (nPartitions < (rank+1) * nMaxLocalPart)
		nLocalPart = std::max(0, static_cast<int>(nPartitions - rank * nMaxLocalPart));

	// Migrate elements and set partition
	logInfo(rank) << "Redistributing elements";

	apf::MeshTag* partitionTag = mesh->createIntTag("partition", 1);

	apf::Migration* plan = new apf::Migration(mesh);

	it = mesh->begin(3);
	for (unsigned int i = 0; i < nLocalElements; i++) {
		apf::MeshEntity* element = mesh->iterate(it);

		int p = part[i];
		mesh->setIntTag(element, partitionTag, &p);

		if (static_cast<unsigned int>(rank) != part[i] / nMaxLocalPart)
			plan->send(element, part[i] / nMaxLocalPart);
	}
	mesh->end(it);

	delete [] part;

	mesh->migrate(plan);

	// Create numbering for the vertices/elements
	apf::GlobalNumbering* vertexNum = apf::makeGlobal(
			apf::numberOwnedNodes(mesh, "vertices"));
	apf::synchronize(vertexNum);
	apf::GlobalNumbering* elementNum = apf::makeGlobal(
			apf::numberOwnedDimension(mesh, "elements", 3));

	// Send/Recv orientation for ghost layers
	logInfo(rank) << "Exchanging ghost layer information";
	PCU_Comm_Begin();
	it = mesh->begin(3);
	while (apf::MeshEntity* element = mesh->iterate(it)) {
		apf::Downward faces;
		mesh->getDownward(element, 2, faces);

		for (unsigned int i = 0; i < 4; i++) {
			if (!mesh->isShared(faces[i]))
				continue;

			// Send local face id and global vertex id
			// of the first vertex
			apf::Copy other = apf::getOtherCopy(mesh, faces[i]);
			PCU_COMM_PACK(other.peer, other.entity);
			int p;
			mesh->getIntTag(element, partitionTag, &p);
			PCU_COMM_PACK(other.peer, p);
			PCU_COMM_PACK(other.peer, i);
			long eid = apf::getNumber(elementNum, apf::Node(element, 0));
			PCU_COMM_PACK(other.peer, eid);
			long vid = getFirstNumOfFace(vertexNum, element, i);
			PCU_COMM_PACK(other.peer, vid);
		}
	}
	mesh->end(it);
	PCU_Comm_Send();

	apf::MeshTag* remoteInfo1 = mesh->createIntTag("remote info 1", 2);
	apf::MeshTag* remoteInfo2 = mesh->createLongTag("remote info 2", 2);

	while (PCU_Comm_Receive()) {
		apf::MeshEntity* face;
		PCU_COMM_UNPACK(face);
		int info1[2];
		PCU_Comm_Unpack(info1, sizeof(info1));
		mesh->setIntTag(face, remoteInfo1, info1);
		long info2[2];
		PCU_Comm_Unpack(info2, sizeof(info2));
		mesh->setLongTag(face, remoteInfo2, info2);
	}

	// Element map is stored as a tag
	apf::MeshTag* localId = mesh->createIntTag("local id", 1);

	// MPI index is also stored as a tag
	apf::MeshTag* mpiIndex = mesh->createIntTag("MPI index", 1);

	// Get tags from mesh sources
	apf::MeshTag* boundaryTag = mesh->findTag("boundary condition");
	apf::MeshTag* groupTag = mesh->findTag("group");

	// Create communicator for I/O
	MPI_Comm commIO;
	MPI_Comm_split(MPI_COMM_WORLD, (nLocalPart > 0 ? 0 : MPI_UNDEFINED), 0, &commIO);

	if (nLocalPart > 0) {
		logInfo(rank) << "Creating element, vertex and boundary maps";
		unsigned int* elemSize = new unsigned int[nLocalPart];
		memset(elemSize, 0, nLocalPart * sizeof(unsigned int));
		std::map<apf::MeshEntity*, unsigned int>* vertexMap =
				new std::map<apf::MeshEntity*, unsigned int>[nLocalPart];
		std::map<unsigned int, std::vector<MPINeighborElement> >* boundaryMap =
				new std::map<unsigned int, std::vector<MPINeighborElement> >[nLocalPart];

		it = mesh->begin(3);
		while (apf::MeshEntity* element = mesh->iterate(it)) {
			// Create element map and count elements
			int p;
			mesh->getIntTag(element, partitionTag, &p);
			int id = elemSize[p % nMaxLocalPart];
			mesh->setIntTag(element, localId, &id);
			elemSize[p % nMaxLocalPart]++;

			// Create vertex map
			apf::Downward vertices;
			mesh->getDownward(element, 0, vertices);
			for (unsigned int i = 0; i < 4; i++) {
				if (vertexMap[p % nMaxLocalPart].find(vertices[i])
						== vertexMap[p % nMaxLocalPart].end()) {
					unsigned int s = vertexMap[p % nMaxLocalPart].size();
					vertexMap[p % nMaxLocalPart][vertices[i]] = s;
				}
			}

			// Compute boundary map
			apf::Downward faces;
			mesh->getDownward(element, 2, faces);

			for (unsigned int i = 0; i < 4; i++) {
				if (mesh->isShared(faces[i])) {
					int info1[2];
					mesh->getIntTag(faces[i], remoteInfo1, info1);
					long info2[2];
					mesh->getLongTag(faces[i], remoteInfo2, info2);
					MPINeighborElement neighbor = {
							element,
							apf::getNumber(elementNum, apf::Node(element, 0)),
							static_cast<int>(i),
							info2[0],
							info1[1]

					};
					boundaryMap[p % nMaxLocalPart][info1[0]].push_back(neighbor);
				} else {
					apf::MeshEntity* n = getFaceElemOppositeElem(mesh, faces[i], element);
					if (!n)
						// Geometric boundary
						continue;

					int np;
					mesh->getIntTag(n, partitionTag, &np);

					if (np == p)
						// Not an MPI boundary
						continue;

					apf::Downward nf;
					mesh->getDownward(n, 2, nf);
					MPINeighborElement neighbor = {
							element,
							apf::getNumber(elementNum, apf::Node(element, 0)),
							static_cast<int>(i),
							apf::getNumber(elementNum, apf::Node(n, 0)),
							apf::findIn(nf, 4, faces[i])

					};
					boundaryMap[p % nMaxLocalPart][np].push_back(neighbor);
				}
			}
		}
		mesh->end(it);

		// Sorting boundary maps and set mpi index tag
		for (unsigned int i = 0; i < nLocalPart; i++) {
			for (std::map<unsigned int, std::vector<MPINeighborElement> >::iterator it = boundaryMap[i].begin();
					it != boundaryMap[i].end(); it++) {
				if (i+(rank*nMaxLocalPart) > it->first)
					std::sort(it->second.begin(), it->second.end(), compareRemoteMPINeighbor);
				else
					std::sort(it->second.begin(), it->second.end(), compareLocalMPINeighbor);

				// After sorting set the mpi index for the elements
				for (int j = 0; j < static_cast<int>(it->second.size()); j++) {
					apf::Downward f;
					mesh->getDownward(it->second[j].element, 2, f);
					mesh->setIntTag(f[it->second[j].localSide], mpiIndex, &j);
				}
			}
		}

		logInfo(rank) << "Computing dimension sizes";
		unsigned int maxSize[4] = {0, 0, 0, 0}; // elements, vertices, boundaries, boundary sizes

		for (unsigned int i = 0; i < nLocalPart; i++) {
			if (elemSize[i] > maxSize[0])
				maxSize[0] = elemSize[i];

			if (vertexMap[i].size() > maxSize[1])
				maxSize[1] = vertexMap[i].size();

			if (boundaryMap[i].size() > maxSize[2])
				maxSize[2] = boundaryMap[i].size();

			for (std::map<unsigned int, std::vector<MPINeighborElement> >::const_iterator it = boundaryMap[i].begin();
					it != boundaryMap[i].end(); it++) {
				if (it->second.size() > maxSize[3])
					maxSize[3] = it->second.size();
			}
		}

		MPI_Allreduce(MPI_IN_PLACE, maxSize, 4, MPI_UNSIGNED, MPI_MAX, commIO);

		logInfo(rank) << "Creating netCDF file";
		// TODO create netCDF file with a new PUML format
		int ncFile;
		checkNcError(nc_create_par(outputFile.c_str(), NC_NETCDF4 | NC_MPIIO,
				commIO, MPI_INFO_NULL, &ncFile));

		// Create netcdf dimensions
		int ncDimDimension;
		nc_def_dim(ncFile, "dimension", 3, &ncDimDimension);

		int ncDimPart;
		nc_def_dim(ncFile, "partitions", nPartitions, &ncDimPart);

		int ncDimElem, ncDimElemSides, ncDimElemVertices;
		nc_def_dim(ncFile, "elements", maxSize[0], &ncDimElem);
		nc_def_dim(ncFile, "element_sides", 4, &ncDimElemSides);
		nc_def_dim(ncFile, "element_vertices", 4, &ncDimElemVertices);

		int ncDimVrtx;
		nc_def_dim(ncFile, "vertices", maxSize[1], &ncDimVrtx);

		int ncDimBnd, ncDimBndElem;
		nc_def_dim(ncFile, "boundaries", maxSize[2], &ncDimBnd);
		nc_def_dim(ncFile, "boundary_elements", maxSize[3], &ncDimBndElem);

		// Create netcdf variables
		int ncVarElemSize;
		checkNcError(nc_def_var(ncFile, "element_size", NC_INT, 1, &ncDimPart, &ncVarElemSize));
		checkNcError(nc_var_par_access(ncFile, ncVarElemSize, NC_COLLECTIVE));

		int ncVarElemVertices;
		int dimsElemVertices[] = {ncDimPart, ncDimElem, ncDimElemVertices};
		checkNcError(nc_def_var(ncFile, "element_vertices", NC_INT, 3, dimsElemVertices, &ncVarElemVertices));
		checkNcError(nc_var_par_access(ncFile, ncVarElemVertices, NC_COLLECTIVE));

		int ncVarElemNeighbors;
		int dimsElemSides[] = {ncDimPart, ncDimElem, ncDimElemSides};
		checkNcError(nc_def_var(ncFile, "element_neighbors", NC_INT, 3, dimsElemSides, &ncVarElemNeighbors));
		checkNcError(nc_var_par_access(ncFile, ncVarElemNeighbors, NC_COLLECTIVE));

		int ncVarElemBoundaries;
		checkNcError(nc_def_var(ncFile, "element_boundaries", NC_INT, 3, dimsElemSides, &ncVarElemBoundaries));
		checkNcError(nc_var_par_access(ncFile, ncVarElemBoundaries, NC_COLLECTIVE));

		int ncVarElemNeighborSides;
		checkNcError(nc_def_var(ncFile, "element_neighbor_sides", NC_INT, 3, dimsElemSides, &ncVarElemNeighborSides));
		checkNcError(nc_var_par_access(ncFile, ncVarElemNeighborSides, NC_COLLECTIVE));

		int ncVarElemSideOrientations;
		checkNcError(nc_def_var(ncFile, "element_side_orientations", NC_INT, 3, dimsElemSides, &ncVarElemSideOrientations));
		checkNcError(nc_var_par_access(ncFile, ncVarElemSideOrientations, NC_COLLECTIVE));

		int ncVarElemNeighborRanks;
		checkNcError(nc_def_var(ncFile, "element_neighbor_ranks", NC_INT, 3, dimsElemSides, &ncVarElemNeighborRanks));
		checkNcError(nc_var_par_access(ncFile, ncVarElemNeighborRanks, NC_COLLECTIVE));

		int ncVarElemMPIIndices;
		checkNcError(nc_def_var(ncFile, "element_mpi_indices", NC_INT, 3, dimsElemSides, &ncVarElemMPIIndices));
		checkNcError(nc_var_par_access(ncFile, ncVarElemMPIIndices, NC_COLLECTIVE));

		int ncVarElemGroup;
		checkNcError(nc_def_var(ncFile, "element_group", NC_INT, 2, dimsElemSides, &ncVarElemGroup));
		checkNcError(nc_var_par_access(ncFile, ncVarElemGroup, NC_COLLECTIVE));

		int ncVarVrtxSize;
		checkNcError(nc_def_var(ncFile, "vertex_size", NC_INT, 1, &ncDimPart, &ncVarVrtxSize));
		checkNcError(nc_var_par_access(ncFile, ncVarVrtxSize, NC_COLLECTIVE));

		int ncVarVrtxCoords;
		int dimsVrtxCoords[] = {ncDimPart, ncDimVrtx, ncDimDimension};
		checkNcError(nc_def_var(ncFile, "vertex_coordinates", NC_DOUBLE, 3, dimsVrtxCoords, &ncVarVrtxCoords));
		checkNcError(nc_var_par_access(ncFile, ncVarVrtxCoords, NC_COLLECTIVE));

		int ncVarBndSize;
		checkNcError(nc_def_var(ncFile, "boundary_size", NC_INT, 1, &ncDimPart, &ncVarBndSize));
		checkNcError(nc_var_par_access(ncFile, ncVarBndSize, NC_COLLECTIVE));

		int ncVarBndElemSize;
		int dimsBndElemSize[] = {ncDimPart, ncDimBnd};
		checkNcError(nc_def_var(ncFile, "boundary_element_size", NC_INT, 2, dimsBndElemSize, &ncVarBndElemSize));
		checkNcError(nc_var_par_access(ncFile, ncVarBndElemSize, NC_COLLECTIVE));

		int ncVarBndElemRank;
		checkNcError(nc_def_var(ncFile, "boundary_element_rank", NC_INT, 2, dimsBndElemSize, &ncVarBndElemRank));
		checkNcError(nc_var_par_access(ncFile, ncVarBndElemRank, NC_COLLECTIVE));

		int ncVarBndElemLocalIds;
		int dimsBndElemLocalIds[] = {ncDimPart, ncDimBnd, ncDimBndElem};
		checkNcError(nc_def_var(ncFile, "boundary_element_localids", NC_INT, 3, dimsBndElemLocalIds, &ncVarBndElemLocalIds));
		checkNcError(nc_var_par_access(ncFile, ncVarBndElemLocalIds, NC_COLLECTIVE));

		checkNcError(nc_enddef(ncFile));

		// Buffers for I/O
		unsigned int* elemVertices = new unsigned int[maxSize[0]*4];
		int* elemNeighbors = new int[maxSize[0]*4];
		int* elemBoundary = new int[maxSize[0]*4];
		unsigned int* elemNbSide = new unsigned int[maxSize[0]*4];
		unsigned int* elemOrientation = new unsigned int[maxSize[0]*4];
		int* elemNbRank = new int[maxSize[0]*4];
		int* elemMPIIndex = new int[maxSize[0]*4];
		int* elemGroup = new int[maxSize[0]];
		double* vertices = new double[maxSize[1]*3];
		int* boundaryIds = new int[maxSize[3]];

		for (unsigned int i = 0; i < nMaxLocalPart; i++) {
			logInfo(rank) << "Writing netCDF file part" << (i+1) << "of" << nMaxLocalPart;

			unsigned int j = i % nLocalPart;

			// Elements
			size_t start[3] = {j + rank*nMaxLocalPart, 0, 0};
			checkNcError(nc_put_var1_uint(ncFile, ncVarElemSize, start, &elemSize[j]));
			size_t count[3] = {1, static_cast<size_t>(elemSize[j]), 4};

			memset(elemBoundary, 0, elemSize[j]*4*sizeof(int));
			memset(elemNbSide, 0, elemSize[j]*4*sizeof(unsigned int));
			memset(elemOrientation, 0, elemSize[j]*4*sizeof(unsigned int));
			memset(elemMPIIndex, 0, elemSize[j]*4*sizeof(int));
			memset(elemGroup, 0, elemSize[j]*sizeof(int));

			it = mesh->begin(3);
			unsigned int k = 0;
			while (apf::MeshEntity* element = mesh->iterate(it)) {
				int part;
				mesh->getIntTag(element, partitionTag, &part);
				if (part % nMaxLocalPart != j)
					continue;

				apf::NewArray<long> vn;
				apf::getElementNumbers(vertexNum, element, vn);

				apf::Downward vert;
				mesh->getDownward(element, 0, vert);
				for (unsigned int l = 0; l < 4; l++)
					elemVertices[k*4 + l] = vertexMap[j][vert[l]];

				apf::Downward faces;
				mesh->getDownward(element, 2, faces);
				for (unsigned int l = 0; l < 4; l++) {
					if (mesh->isShared(faces[l])) {
						// Partition boundary on another process
						elemNeighbors[k*4 + FACE2INTERNAL[l]] = elemSize[j];

						int info1[2];
						mesh->getIntTag(faces[l], remoteInfo1, info1);
						elemNbSide[k*4 + FACE2INTERNAL[l]] = FACE2INTERNAL[info1[1]];

						long info2[2];
						mesh->getLongTag(faces[l], remoteInfo2, info2);
						elemOrientation[k*4 + FACE2INTERNAL[l]] = getOrientation(&vn[0], l, info2[1]);

						elemNbRank[k*4 + FACE2INTERNAL[l]] = info1[0];

						mesh->getIntTag(faces[l], mpiIndex, &elemMPIIndex[k*4 + FACE2INTERNAL[l]]);
					} else {
						apf::MeshEntity* n = getFaceElemOppositeElem(mesh, faces[l], element);

						if (n == 0L) {
							// Geometric boundary
							elemNeighbors[k*4 + FACE2INTERNAL[l]] = elemSize[j];

							elemNbRank[k*4 + FACE2INTERNAL[l]] = part;
						} else {
							int np;
							mesh->getIntTag(n, partitionTag, &np);

							if (np % nMaxLocalPart == j) {
								// Same partition
								mesh->getIntTag(n, localId, &elemNeighbors[k*4 + FACE2INTERNAL[l]]);
							} else {
								// Partition boundary on same process
								elemNeighbors[k*4 + FACE2INTERNAL[l]] = elemSize[j];

								mesh->getIntTag(faces[l], mpiIndex, &elemMPIIndex[k*4 + FACE2INTERNAL[l]]);
							}

							apf::Downward nf;
							mesh->getDownward(n, 2, nf);
							int nfId = apf::findIn(nf, 4, faces[l]);
							elemNbSide[k*4 + FACE2INTERNAL[l]] = FACE2INTERNAL[nfId];

							elemOrientation[k*4 + FACE2INTERNAL[l]] =
									getOrientation(&vn[0], l, getFirstNumOfFace(vertexNum, n, nfId));

							mesh->getIntTag(n, partitionTag, &elemNbRank[k*4 + FACE2INTERNAL[l]]);
						}
					}

					if (boundaryTag && mesh->hasTag(faces[l], boundaryTag))
						mesh->getIntTag(faces[l], boundaryTag, &elemBoundary[k*4 + FACE2INTERNAL[l]]);

				}

				if (groupTag && mesh->hasTag(element, groupTag))
					mesh->getIntTag(element, groupTag, &elemGroup[k]);

				k++;
			}

			checkNcError(nc_put_vara_uint(ncFile, ncVarElemVertices, start, count, elemVertices));
			checkNcError(nc_put_vara_int(ncFile, ncVarElemNeighbors, start, count, elemNeighbors));
			checkNcError(nc_put_vara_int(ncFile, ncVarElemBoundaries, start, count, elemBoundary));
			checkNcError(nc_put_vara_uint(ncFile, ncVarElemNeighborSides, start, count, elemNbSide));
			checkNcError(nc_put_vara_uint(ncFile, ncVarElemSideOrientations, start, count, elemOrientation));
			checkNcError(nc_put_vara_int(ncFile, ncVarElemNeighborRanks, start, count, elemNbRank));
			checkNcError(nc_put_vara_int(ncFile, ncVarElemMPIIndices, start, count, elemMPIIndex));
			checkNcError(nc_put_vara_int(ncFile, ncVarElemGroup, start, count, elemGroup));

			// Vertices
			unsigned int vertexSize = vertexMap[j].size();
			checkNcError(nc_put_var1_uint(ncFile, ncVarVrtxSize, start, &vertexSize));
			count[1] = vertexSize; count[2] = 3;
			for (std::map<apf::MeshEntity*, unsigned int>::const_iterator it = vertexMap[j].begin();
					it != vertexMap[j].end(); it++) {
				apf::Vector3 point;
				mesh->getPoint(it->first, 0, point);
				point.toArray(&vertices[it->second*3]);
			}
			checkNcError(nc_put_vara_double(ncFile, ncVarVrtxCoords, start, count, vertices));

			// Boundaries
			int s = boundaryMap[j].size();
			checkNcError(nc_put_var1_int(ncFile, ncVarBndSize, start, &s));

			unsigned int bndCount = 0;
			int size;	// need to declare them outside of the loop (to reuse last values during
			int remoteRank;	// collective I/O)

			count[0] = 1; count[1] = 1;
			for (std::map<unsigned int, std::vector<MPINeighborElement> >::const_iterator it = boundaryMap[j].begin();
					it != boundaryMap[j].end(); it++) {
				start[1] = bndCount;

				size = it->second.size();
				checkNcError(nc_put_var1_int(ncFile, ncVarBndElemSize, start, &size));
				remoteRank = it->first;
				checkNcError(nc_put_var1_int(ncFile, ncVarBndElemRank, start, &remoteRank));

				// Fill local ids
				for (size_t k = 0; k < it->second.size(); k++)
					mesh->getIntTag(it->second[k].element, localId, &boundaryIds[k]);

				count[2] = size;
				checkNcError(nc_put_vara_int(ncFile, ncVarBndElemLocalIds, start, count, boundaryIds));

				bndCount++;
			}

			// For collective I/O
			for (; bndCount < maxSize[3]; bndCount++) {
				checkNcError(nc_put_var1_int(ncFile, ncVarBndElemSize, start, &size));
				checkNcError(nc_put_var1_int(ncFile, ncVarBndElemRank, start, &remoteRank));

				checkNcError(nc_put_vara_int(ncFile, ncVarBndElemLocalIds, start, count, boundaryIds));
			}
		}

		delete [] elemVertices;
		delete [] elemNeighbors;
		delete [] elemBoundary;
		delete [] elemNbSide;
		delete [] elemNbRank;
		delete [] elemMPIIndex;
		delete [] vertices;
		delete [] boundaryIds;

		delete [] elemSize;
		delete [] vertexMap;
		delete [] boundaryMap;

		checkNcError(nc_close(ncFile));

		MPI_Comm_free(&commIO);
	}

	delete meshInput;

	logInfo(rank) << "Finished successfully";

	PCU_Comm_Free();

	MPI_Finalize();
	return 0;
}
