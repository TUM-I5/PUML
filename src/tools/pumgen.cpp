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

#include "utils/args.h"
#include "utils/logger.h"

#include <parmetis.h>

#include <netcdf_par.h>
#include <netcdf.h>

#include "input/Apf.h"
#include "input/SerialMeshFile.h"
#include "input/ApfNative.h"
#ifdef USE_SIMMOD
#include "input/SimModSuite.h"
#endif // USE_SIMMOD
#include "meshreader/ParallelGambitReader.h"

struct Element {
	/** Global id */
	unsigned int id;
	/** Vertex ids */
	unsigned int vertex[4];
	/** The group number */
	unsigned int group;
	/** The boundary conditions */
	unsigned int boundaries[4];
	/** The neighbor ids */
	int neighbors[4];

	/**
	 * Compare elements by the global id
	 */
	bool operator<(const Element& e) const
	{
		return id < e.id;
	}
};

struct MPINeighborElement {
	/** Local number of the local element */
	int localElement;
	/** Side of the local element */
	int localSide;
	/** Global number neighbor element */
	int neighborElement;
	/** Side of the neighbor element */
	int neighborSide;
};

void checkNcError(int error)
{
	if (error != NC_NOERR)
		logError() << "Error while writing netCDF file:" << nc_strerror(error);
}

bool compareLocalMPINeighbor(const MPINeighborElement &elem1, const MPINeighborElement &elem2)
{
	return (elem1.localElement < elem2.localElement)
			|| (elem1.localElement == elem2.localElement && elem1.localSide < elem2.localSide);
}

bool compareRemoteMPINeighbor(const MPINeighborElement &elem1, const MPINeighborElement &elem2)
{
	return (elem1.neighborElement < elem2.neighborElement)
				|| (elem1.neighborElement == elem2.neighborElement && elem1.neighborSide < elem2.neighborSide);
}

/**
 * @param[out] result Upon return a bit mask that indicates which values where found and which not
 */
template<class InputIterator1, class InputIterator2>
void findAllBit(InputIterator1 first, InputIterator1 last, InputIterator2 firstVal, InputIterator2 lastVal, int* result)
{
	size_t pos = 0;
	while (firstVal != lastVal) {
		if (pos % (sizeof(int)*8) == 0)
			result[pos / (sizeof(int)*8)] = 0;

		InputIterator1 ret = std::find(first, last, *firstVal);
		if (ret != last)
			result[pos / (sizeof(int)*8)] |= 1<<(pos % (sizeof(int)*8));

		firstVal++;
		pos++;
	}
}

/**
 * @return The corresponding face to the bitmask
 */
int faceOfMask(int mask)
{
	switch (mask) {
	case 7:
		return 0;
	case 11:
		return 1;
		break;
	case 13:
		return 2;
	case 14:
		return 3;
	default:
		logError() << "Dual graph is wrong";
	}

	return -1;
}

/** Maps from element + face to orientation */
const static int face2orientation[4][4] = {
		{0, 0, 0, -1},
		{2, 1, -1, 0},
		{1, -1, 2, 1},
		{-1, 2, 1, 2}
		//{0, 2, 1},
		//{0, 1, 3},
		//{0, 3, 2},
		//{1, 2, 3}
};

const static int METIS_RANDOM_SEED = 42;

int main(int argc, char* argv[])
{
	int rank = 0;
	int processes = 1;

	int threadSupport;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &threadSupport);

	if (threadSupport != MPI_THREAD_SINGLE)
		// I'm not sure yet if we really need thread support
		logError() << "Threading is not supported in this MPI version";

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &processes);

	PCU_Comm_Init();

	// Parse command line arguments
	utils::Args args;
	const char* source[] = {"gambit", "apf", "simmodsuite"};
	args.addEnumOption("source", source, 's', "Mesh source (default: gambit)", false);
	args.addOption("dump", 'd', "Dump the mesh before partitioning it",
			utils::Args::Required, false);
	args.addOption("model", 0, "Dump/Load a specific model file",
			utils::Args::Required, false);
	args.addOption("license", 'l', "License file (only used by SimModSuite)",
			utils::Args::Required, false);
	args.addOption("cad", 'c', "CAD file (only used by SimModSuite)",
			utils::Args::Required, false);
	args.addOption("mesh", 0, "Mesh attributes name (only used by SimModSuite, default: \"mesh\")",
			utils::Args::Required, false);
	args.addOption("analysis", 0, "Analysis attributes name (only used by SimModSuite, default: \"analysis\")",
			utils::Args::Required, false);
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

	if (nPartitions <= ((nPartitions + processes - 1) / processes) * (processes-1))
		logError() << "Not every process will get at least one partition, use a smaller number of ranks";

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
	switch (args.getArgument<int>("source", 0)) {
	case 0:
		logInfo(rank) << "Using Gambit mesh";
		meshInput = new SerialMeshFile<ParallelGambitReader>(inputFile);
		break;
	case 1:
		logInfo(rank) << "Using APF native format";
		meshInput = new ApfNative(inputFile);
		break;
	case 2:
#ifdef USE_SIMMOD
		logInfo(rank) << "Using SimModSuite";

		meshInput = new SimModSuite(inputFile,
				args.getArgument<const char*>("cad", 0L),
				args.getArgument<const char*>("license", 0L),
				args.getArgument<const char*>("mesh", "mesh"),
				args.getArgument<const char*>("analysis", "analysis"),
				args.getArgument<int>("enforce-size", 0));

#else // USE_SIMMOD
		logError() << "SimModSuite is not supported in this version";
#endif // USE_SIMMOD
		break;
	default:
		logError() << "Unknown source";
	}

	Apf* mesh = new Apf(meshInput->getMesh());

	const char* dumpFile = args.getArgument<const char*>("dump", 0L);
	if (dumpFile)
		mesh->write(dumpFile, args.getArgument<const char*>("model", 0L));

	delete meshInput;

	// Read elements
	// TODO only tetrahedral meshes are currently supported
	unsigned int* elements;
	MPI_Alloc_mem(mesh->nLocalElements()*4*sizeof(unsigned int), MPI_INFO_NULL, &elements);
	mesh->getElements(elements);

	// TODO wrap this code part so we can switch to different libraries
	logInfo(rank) << "Creating dual graph with METIS";
	idx_t* elemDist = new idx_t[processes+1];
	for (int i = 0; i < processes; i++)
		elemDist[i] = mesh->elemStart(i);
	elemDist[processes] = mesh->nElements();

	idx_t* eptr = new idx_t[mesh->nLocalElements()+1];
	for (unsigned int i = 0; i < mesh->nLocalElements()+1; i++)
		eptr[i] = i * 4;

	idx_t* eind = new idx_t[mesh->nLocalElements()*4];
	for (unsigned int i = 0; i < mesh->nLocalElements()*4; i++)
		eind[i] = elements[i];

	idx_t numflag = 0;
	idx_t ncommonnodes = 3;

	idx_t* xadj;
	idx_t* adjncy;
	MPI_Comm commWorld = MPI_COMM_WORLD;
	if (ParMETIS_V3_Mesh2Dual(elemDist, eptr, eind, &numflag, &ncommonnodes, &xadj, &adjncy, &commWorld) != METIS_OK)
		logError() << "Could not create dual graph";

	delete [] eptr;
	delete [] eind;

	logInfo(rank) << "Creating neighborhood information";
	int* elementNeighbors = new int[mesh->nLocalElements()*4];
	for (unsigned int i = 0; i < mesh->nLocalElements(); i++) {
		idx_t j;
		for (j = 0; j < xadj[i+1]-xadj[i]; j++)
			elementNeighbors[i*4 + j] = adjncy[xadj[i] + j];

		for (; j < 4; j++)
			elementNeighbors[i*4 + j] = -1;

	}

	// Read boundary information from the mesh
	unsigned int* boundaries = new unsigned int[mesh->nLocalElements()*4];
	memset(boundaries, 0, mesh->nLocalElements()*4*sizeof(unsigned int));
	mesh->getBoundaries(boundaries);

	logInfo(rank) << "Creating partitions with METIS";
	idx_t wgtflag = 0;
	idx_t ncon = 1;
	idx_t nparts = nPartitions;

	real_t* tpwgts = new real_t[nPartitions];
	for (unsigned int i = 0; i < nPartitions; i++)
		tpwgts[i] = 1./nPartitions;

	real_t ubev = 1.01;
	idx_t options[3] = {1, 0, METIS_RANDOM_SEED};

	idx_t edgecut;
	idx_t* part = new idx_t[mesh->nLocalElements()];

	if (ParMETIS_V3_PartKway(elemDist, xadj, adjncy, 0L, 0L, &wgtflag, &numflag, &ncon,
			&nparts, tpwgts, &ubev, options, &edgecut, part, &commWorld) != METIS_OK)
		logError() << "Could create partitions";

#if 0
	// Test partitioning
	std::stringstream ss;
	ss << "test-" << rank << ".data";

	std::ofstream output(ss.str());

	for (unsigned int i = 0; i < nLocalElements; i++) {
		output << (i+rank*nMaxLocalElements) << ' ' << part[i] << std::endl;
	}
#endif

	delete [] elemDist;
	delete [] tpwgts;

	METIS_Free(xadj);
	METIS_Free(adjncy);

	// Read group information
	unsigned int* elementGroups = new unsigned int[mesh->nLocalElements()];
	mesh->getGroups(elementGroups);

	logInfo(rank) << "Redistributing elements";
	unsigned int nMaxLocalPart = (nPartitions + processes - 1) / processes;
	unsigned int nLocalPart = nMaxLocalPart;
	if (rank == processes - 1)
		nLocalPart = nPartitions - (processes-1) * nMaxLocalPart;

	// Redistribute the elements among the processors according to there partitions
	// Count number local elements in each partition
	unsigned int* localPartPtr = new unsigned int[nLocalPart+1];
	localPartPtr[0] = 0;

	Element* localPartElements;

	unsigned int pos;	// Helper variable

	std::vector<unsigned int>* elementPerPart = new std::vector<unsigned int>[nPartitions];
	for (unsigned int i = 0; i < mesh->nLocalElements(); i++)
		elementPerPart[part[i]].push_back(i);
	unsigned int* nElementPerPart = new unsigned int[nPartitions];
	for (unsigned int i = 0; i < nPartitions; i++)
		nElementPerPart[i] = elementPerPart[i].size();

	// Tell each processor how many elements we are going to send
	int* sendcounts = new int[processes];
	int* senddispls = new int[processes];
	for (int i = 0; i < processes-1; i++) {
		sendcounts[i] = nMaxLocalPart;
		senddispls[i] = nMaxLocalPart * i;
	}
	sendcounts[processes-1] = nPartitions - (processes-1) * nMaxLocalPart;
	senddispls[processes-1] = (processes-1) * nMaxLocalPart;

	unsigned int* recvSize = new unsigned int[processes * nLocalPart];
	int* recvcounts = new int[processes];
	int* recvdispls = new int[processes];
	for (int i = 0; i < processes; i++) {
		recvcounts[i] = nLocalPart;
		recvdispls[i] = nLocalPart * i;
	}

	MPI_Alltoallv(nElementPerPart, sendcounts, senddispls, MPI_UNSIGNED,
			recvSize, recvcounts, recvdispls, MPI_UNSIGNED, MPI_COMM_WORLD);

	// Now send the the element to all processes
	Element* sendElements = new Element[mesh->nLocalElements()];
	pos = 0;
	for (unsigned int i = 0; i < nPartitions; i++) {
		for (unsigned int j = 0; j < nElementPerPart[i]; j++) {
			// Set global id
			sendElements[pos].id = elementPerPart[i][j] + mesh->elemStart(rank);
			memcpy(sendElements[pos].vertex,
					&elements[elementPerPart[i][j]*4], 4*sizeof(unsigned int));
			sendElements[pos].group = elementGroups[elementPerPart[i][j]];
			memcpy(sendElements[pos].boundaries,
					&boundaries[elementPerPart[i][j]*4], 4*sizeof(unsigned int));
			memcpy(sendElements[pos].neighbors,
					&elementNeighbors[elementPerPart[i][j]*4], 4*sizeof(int));
			pos++;
		}
	}

	delete [] elementGroups;
	delete [] boundaries;
	delete [] elementNeighbors;

	memset(sendcounts, 0, sizeof(int) * processes);
	for (unsigned int i = 0; i < nPartitions; i++) {
		sendcounts[i / nMaxLocalPart] += nElementPerPart[i];
	}
	senddispls[0] = 0;
	for (int i = 1; i < processes; i++)
		senddispls[i] = senddispls[i-1] + sendcounts[i-1];

	memset(recvcounts, 0, sizeof(int) * processes);
	for (int i = 0; i < processes; i++) {
		for (unsigned int j = 0; j < nLocalPart; j++)
			recvcounts[i] += recvSize[i*nLocalPart + j];
	}
	recvdispls[0] = 0;
	for (int i = 1; i < processes; i++)
		recvdispls[i] = recvdispls[i-1] + recvcounts[i-1];

	unsigned int recvSum = 0;
	for (int i = 0; i < processes; i++) {
		recvSum += recvcounts[i];
	}
	Element* recvElements = new Element[recvSum];

	// Create the MPI datatype
	int blockLength[] = {1, 4, 1, 4, 4};
	MPI_Aint displacement[] = {0, reinterpret_cast<uintptr_t>(recvElements[0].vertex)
			- reinterpret_cast<uintptr_t>(&recvElements[0]),
			reinterpret_cast<uintptr_t>(&recvElements[0].group)
			- reinterpret_cast<uintptr_t>(&recvElements[0]),
			reinterpret_cast<uintptr_t>(recvElements[0].boundaries)
			- reinterpret_cast<uintptr_t>(&recvElements[0]),
			reinterpret_cast<intptr_t>(recvElements[0].neighbors)
			- reinterpret_cast<intptr_t>(&recvElements[0])};
	MPI_Datatype type[] = {MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED,
			MPI_UNSIGNED, MPI_INT};
	MPI_Datatype elementType;
	MPI_Type_create_struct(5, blockLength, displacement, type, &elementType);
	MPI_Type_commit(&elementType);

	MPI_Alltoallv(sendElements, sendcounts, senddispls, elementType,
			recvElements, recvcounts, recvdispls, elementType, MPI_COMM_WORLD);

	MPI_Type_free(&elementType);

	delete [] sendcounts;
	delete [] senddispls;
	delete [] recvcounts;
	delete [] recvdispls;
	delete [] sendElements;

	delete [] nElementPerPart;
	delete [] elementPerPart;

	// Compute the size of each local partition
	unsigned int* partSize = new unsigned int[nLocalPart];
	memset(partSize, 0, sizeof(unsigned int)*nLocalPart);
	for (int i = 0; i < processes; i++) {
		for (unsigned int j = 0; j < nLocalPart; j++) {
			partSize[j] += recvSize[i*nLocalPart + j];
		}
	}
	for (unsigned int i = 1; i < nLocalPart+1; i++) {
		localPartPtr[i] = localPartPtr[i-1] + partSize[i-1];
	}

	delete [] partSize;

	// Sort elements by partitions
	localPartElements = new Element[localPartPtr[nLocalPart]];
	unsigned int* actPartPos = new unsigned int[nLocalPart];
	memcpy(actPartPos, localPartPtr, sizeof(unsigned int)*nLocalPart);
	unsigned int actRecvPos = 0;
	for (int i = 0; i < processes; i++) {
		for (unsigned int j = 0; j < nLocalPart; j++) {
			memcpy(&localPartElements[actPartPos[j]], &recvElements[actRecvPos],
					sizeof(Element)*recvSize[i*nLocalPart + j]);
			actPartPos[j] += recvSize[i*nLocalPart + j];
			actRecvPos += recvSize[i*nLocalPart + j];
		}
	}
	delete [] actPartPos;

	delete [] recvElements;
	delete [] recvSize;

	// Save the partition information because other processes still need it
	// But copy it to an unsigned int
	unsigned int* elementPart;
	MPI_Alloc_mem(mesh->nLocalElements()*sizeof(unsigned int), MPI_INFO_NULL, &elementPart);
	for (unsigned int i = 0; i < mesh->nLocalElements(); i++)
		elementPart[i] = part[i];

	delete [] part;

	logInfo(rank) << "Sorting element ids";
	for (unsigned int i = 0; i < nLocalPart; i++)
		std::sort(&localPartElements[localPartPtr[i]],
				&localPartElements[localPartPtr[i]]+localPartPtr[i+1]-localPartPtr[i]);

	logInfo(rank) << "Computing element data";
	// Create MPI window to allow access to all elements
	MPI_Win elementWindow;
	MPI_Win_create(elements, mesh->nLocalElements()*4*sizeof(unsigned int),
			sizeof(unsigned int), MPI_INFO_NULL, MPI_COMM_WORLD, &elementWindow);

	// Create window for partition information
	MPI_Win elementPartWindow;
	MPI_Win_create(elementPart, mesh->nLocalElements()*sizeof(unsigned int),
			sizeof(unsigned int), MPI_INFO_NULL, MPI_COMM_WORLD, &elementPartWindow);

	// Create map from local elements to ids inside their partition
	std::map<unsigned int, unsigned int> globalElem2PartId;

	unsigned int curPart = 0;
	for (unsigned int i = 0; i < localPartPtr[nLocalPart]; i++) {
		if (i >= localPartPtr[curPart+1])
			curPart++;

		globalElem2PartId[localPartElements[i].id] = i-localPartPtr[curPart];
	}

	int* localPartElemNb = new int[localPartPtr[nLocalPart]*4];
	for (unsigned int i = 0; i < localPartPtr[nLocalPart]*4; i++)
		localPartElemNb[i] = -1;
	unsigned int* localPartElemNbSide = new unsigned int[localPartPtr[nLocalPart]*4];
	memset(localPartElemNbSide, 0, sizeof(unsigned int)*4*localPartPtr[nLocalPart]);
	unsigned int* localPartElemNbOrient = new unsigned int[localPartPtr[nLocalPart]*4];
	memset(localPartElemNbOrient, 0, sizeof(unsigned int)*4*localPartPtr[nLocalPart]);
	unsigned int* localPartElemBoundary = new unsigned int[localPartPtr[nLocalPart]*4];
	int* localPartElemRanks = new int[localPartPtr[nLocalPart]*4];
	std::map<unsigned int, std::vector<MPINeighborElement> >* boundaryMaps
		= new std::map<unsigned int, std::vector<MPINeighborElement> >[nLocalPart];
	for (unsigned int i = 0; i < localPartPtr[nLocalPart]*4; i++)
		localPartElemRanks[i] = -1;
	curPart = 0;
	for (unsigned int i = 0; i < localPartPtr[nLocalPart]; i++) {
		if (i >= localPartPtr[curPart+1])
			curPart++;

		// Set boundary (do this for all boundaries, because of inner boundary conditions)
		memcpy(&localPartElemBoundary[i*4], localPartElements[i].boundaries, 4*sizeof(unsigned int));

		for (int j = 0; j < 4; j++) {
			if (localPartElements[i].neighbors[j] < 0)
				continue;

			// Check whether we are in the same partition
			unsigned int neighborPart;
			MPI_Win_lock(MPI_LOCK_SHARED, mesh->rankOfElem(localPartElements[i].neighbors[j]),
					MPI_MODE_NOCHECK, elementPartWindow);
			MPI_Get(&neighborPart, 1, MPI_UNSIGNED,
					mesh->rankOfElem(localPartElements[i].neighbors[j]),
					mesh->posOfElem(localPartElements[i].neighbors[j]),
					1, MPI_UNSIGNED, elementPartWindow);
			MPI_Win_unlock(mesh->rankOfElem(localPartElements[i].neighbors[j]), elementPartWindow);

			unsigned int neighborVertices[4];
			MPI_Win_lock(MPI_LOCK_SHARED, mesh->rankOfElem(localPartElements[i].neighbors[j]),
					MPI_MODE_NOCHECK, elementWindow);
			MPI_Get(neighborVertices, 4, MPI_UNSIGNED,
					mesh->rankOfElem(localPartElements[i].neighbors[j]),
					mesh->posOfElem(localPartElements[i].neighbors[j]) * 4,
					4, MPI_UNSIGNED, elementWindow);
			MPI_Win_unlock(mesh->rankOfElem(localPartElements[i].neighbors[j]), elementWindow);

			int mask;
			findAllBit(neighborVertices, neighborVertices+4, localPartElements[i].vertex, localPartElements[i].vertex+4, &mask);
			int face = faceOfMask(mask);

			// Set neighbor partition
			localPartElemRanks[i*4 + face] = neighborPart;

			// Find the side of the neighbor element
			findAllBit(localPartElements[i].vertex, localPartElements[i].vertex+4, neighborVertices, neighborVertices+4, &mask);
			int face2 = faceOfMask(mask);
			localPartElemNbSide[i*4 + face] = face2;

			if (neighborPart == curPart+(rank*nMaxLocalPart))
				// Set neighbor only if its the same partition
				localPartElemNb[i*4 + face] = globalElem2PartId.at(localPartElements[i].neighbors[j]);
			else {
				// Mark id as neighbor
				MPINeighborElement neighborElement = {static_cast<int>(i-localPartPtr[curPart]),
						face, localPartElements[i].neighbors[j], face2};
				boundaryMaps[curPart][neighborPart].push_back(neighborElement);
			}

			// Find the first vertex of the local face in the neighboring face
			int orientation = face2orientation[std::find(neighborVertices, neighborVertices+4,
					localPartElements[i].vertex[face < 3 ? 0 : 1]) - neighborVertices][face2];
			if (orientation < 0)
				logError() << "Dual mesh is incorrect";
			localPartElemNbOrient[i*4 + face] = orientation;
		}

		// Set all not local neighbors
		for (int j = 0; j < 4; j++) {
			if (localPartElemNb[i*4 + j] < 0) {
				localPartElemNb[i*4 + j] = localPartPtr[curPart+1]-localPartPtr[curPart];

				// Set rank on real boundaries
				if (localPartElemRanks[i*4 + j] < 0)
					localPartElemRanks[i*4 + j] = curPart+(rank*nMaxLocalPart);
			}
		}
	}

	logInfo(rank) << "Sorting MPI boundaries";
	unsigned int* localPartElemMPI = new unsigned int[localPartPtr[nLocalPart]*4];
	memset(localPartElemMPI, 0, sizeof(unsigned int)*localPartPtr[nLocalPart]*4);
	#pragma omp parallel for
	for (unsigned int i = 0; i < nLocalPart; i++) {
		for (std::map<unsigned int, std::vector<MPINeighborElement> >::iterator j = boundaryMaps[i].begin();
				j != boundaryMaps[i].end(); j++) {
			if (i+(rank*nMaxLocalPart) > j->first)
				std::sort(j->second.begin(), j->second.end(), compareRemoteMPINeighbor);
			else
				std::sort(j->second.begin(), j->second.end(), compareLocalMPINeighbor);

			// After sorting set the mpi index for the elements
			for (unsigned int k = 0; k < j->second.size(); k++)
				localPartElemMPI[(j->second[k].localElement + localPartPtr[i])*4 + j->second[k].localSide] = k;
		}
	}

	logInfo(rank) << "Creating vertices maps";
	std::map<unsigned int, unsigned int>* localPartVerticesMaps = new std::map<unsigned int, unsigned int>[nLocalPart];
	unsigned int* localPartVertexSize = new unsigned int[nLocalPart];
	memset(localPartVertexSize, 0, sizeof(unsigned int)*nLocalPart);
	#pragma omp parallel for
	for (unsigned int i = 0; i < nLocalPart; i++) {
		for (unsigned int j = localPartPtr[i]; j < localPartPtr[i+1]; j++) {
			for (int k = 0; k < 4; k++) {
				if (localPartVerticesMaps[i].find(localPartElements[j].vertex[k]) == localPartVerticesMaps[i].end()) {
					// New vertices
					localPartVerticesMaps[i][localPartElements[j].vertex[k]] = localPartVertexSize[i];
					localPartElements[j].vertex[k] = localPartVertexSize[i];
					localPartVertexSize[i]++;
				} else {
					// Vertex already found
					localPartElements[j].vertex[k] = localPartVerticesMaps[i][localPartElements[j].vertex[k]];
				}
			}
		}
	}

	logInfo(rank) << "Computing maximum partition size";
	// Compute max number of elements in a partition
	unsigned int maxPartSize = localPartPtr[1];
	for (unsigned int i = 1; i < nLocalPart; i++)
		if (maxPartSize < localPartPtr[i+1]-localPartPtr[i])
			maxPartSize = localPartPtr[i+1]-localPartPtr[i];
	MPI_Allreduce(MPI_IN_PLACE, &maxPartSize, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);

	// Compute max number of vertices in a partition
	unsigned int maxVertices = localPartVertexSize[0];
	for (unsigned int i = 1; i < nLocalPart; i++)
		if (maxVertices < localPartVertexSize[i])
			maxVertices = localPartVertexSize[i];
	MPI_Allreduce(MPI_IN_PLACE, &maxVertices, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);

	// Compute max number of boundaries per partition
	unsigned int maxBoundaries = boundaryMaps[0].size();
	for (unsigned int i = 1; i < nLocalPart; i++)
		if (maxBoundaries < boundaryMaps[i].size())
			maxBoundaries = boundaryMaps[i].size();
	MPI_Allreduce(MPI_IN_PLACE, &maxBoundaries, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);

	// Compute max boundary size
	unsigned int maxBoundarySize = 0;
	for (unsigned int i = 0; i < nLocalPart; i++) {
		for (std::map<unsigned int, std::vector<MPINeighborElement> >::const_iterator j = boundaryMaps[i].begin();
				j != boundaryMaps[i].end(); j++)
			if (maxBoundarySize < j->second.size())
				maxBoundarySize = j->second.size();

	}
	MPI_Allreduce(MPI_IN_PLACE, &maxBoundarySize, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);

	// Create element vertex buffer
	unsigned int* localPartElemVrtx = new unsigned int[localPartPtr[nLocalPart]*4];
	for (unsigned int i = 0; i < localPartPtr[nLocalPart]; i++)
		memcpy(&localPartElemVrtx[i*4], localPartElements[i].vertex, sizeof(unsigned int)*4);

	// Create element group buffer
	unsigned int* localPartElemGroup = new unsigned int[localPartPtr[nLocalPart]];
	for (unsigned int i = 0; i < localPartPtr[nLocalPart]; i++)
		localPartElemGroup[i] = localPartElements[i].group;

	logInfo(rank) << "Creating netCDF file";
	// TODO create netCDF file with a new PUML format
	int ncFile;
	checkNcError(nc_create_par(outputFile.c_str(), NC_NETCDF4 | NC_MPIIO, MPI_COMM_WORLD, MPI_INFO_NULL, &ncFile));

	// Create netcdf dimensions
	int ncDimDimension;
	nc_def_dim(ncFile, "dimension", 3, &ncDimDimension);

	int ncDimPart;
	nc_def_dim(ncFile, "partitions", nPartitions, &ncDimPart);

	int ncDimElem, ncDimElemSides, ncDimElemVertices;
	nc_def_dim(ncFile, "elements", maxPartSize, &ncDimElem);
	nc_def_dim(ncFile, "element_sides", 4, &ncDimElemSides);
	nc_def_dim(ncFile, "element_vertices", 4, &ncDimElemVertices);

	int ncDimVrtx;
	nc_def_dim(ncFile, "vertices", maxVertices, &ncDimVrtx);

	int ncDimBnd, ncDimBndElem;
	nc_def_dim(ncFile, "boundaries", maxBoundaries, &ncDimBnd);
	nc_def_dim(ncFile, "boundary_elements", maxBoundarySize, &ncDimBndElem);

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

	logInfo(rank) << "Writing element data";
	for (unsigned int i = 0; i < nMaxLocalPart; i++) {
		unsigned int j = i % nLocalPart;

		size_t start[3] = {j + rank*nMaxLocalPart, 0, 0};
		int size = localPartPtr[j+1]-localPartPtr[j];
		checkNcError(nc_put_var1_int(ncFile, ncVarElemSize, start, &size));
		size_t count[3] = {1, static_cast<size_t>(size), 4};
		checkNcError(nc_put_vara_uint(ncFile, ncVarElemVertices, start, count, &localPartElemVrtx[localPartPtr[j]*4]));
		checkNcError(nc_put_vara_int(ncFile, ncVarElemNeighbors, start, count, &localPartElemNb[localPartPtr[j]*4]));
		checkNcError(nc_put_vara_uint(ncFile, ncVarElemBoundaries, start, count, &localPartElemBoundary[localPartPtr[j]*4]));
		checkNcError(nc_put_vara_uint(ncFile, ncVarElemNeighborSides, start, count, &localPartElemNbSide[localPartPtr[j]*4]));
		checkNcError(nc_put_vara_uint(ncFile, ncVarElemSideOrientations, start, count, &localPartElemNbOrient[localPartPtr[j]*4]));
		checkNcError(nc_put_vara_int(ncFile, ncVarElemNeighborRanks, start, count, &localPartElemRanks[localPartPtr[j]*4]));
		checkNcError(nc_put_vara_uint(ncFile, ncVarElemMPIIndices, start, count, &localPartElemMPI[localPartPtr[j]*4]));
		checkNcError(nc_put_vara_uint(ncFile, ncVarElemGroup, start, count, &localPartElemGroup[localPartPtr[j]]));
	}

	// Delete data we have already written to disk
	delete [] localPartElemVrtx;
	delete [] localPartElemNb;
	delete [] localPartElemNbSide;
	delete [] localPartElemNbOrient;
	delete [] localPartElemBoundary;
	delete [] localPartElemRanks;
	delete [] localPartElemMPI;
	delete [] localPartElemGroup;

	MPI_Win_free(&elementWindow);
	MPI_Win_free(&elementPartWindow);
	MPI_Free_mem(elements);
	MPI_Free_mem(elementPart);
	delete [] localPartPtr;
	delete [] localPartElements;

	// Load vertices
	double* vertices;
	MPI_Alloc_mem(mesh->nLocalVertices()*3*sizeof(double), MPI_INFO_NULL, &vertices);
	mesh->getVertices(vertices);

	// Create MPI window to get access to all vertices
	MPI_Win verticesWindow;
	MPI_Win_create(vertices, mesh->nLocalVertices()*3*sizeof(double),
			sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &verticesWindow);

	logInfo(rank) << "Writing vertices data";
	double* localVertices = new double[maxVertices*3];
	bool* localVerticesTransfered = new bool[maxVertices]; // True if we already have the corresponding vertex

	for (unsigned int i = 0; i < nMaxLocalPart; i++) {
		if (i < nLocalPart) {
			// Only get the coordinates if this is a real partition
			memset(localVerticesTransfered, 0, sizeof(bool)*maxVertices);

			for (std::map<unsigned int, unsigned int>::const_iterator j = localPartVerticesMaps[i].begin();
					j != localPartVerticesMaps[i].end(); j++)  {
				// We already have this vertex?
				if (localVerticesTransfered[j->second])
					continue;

				// Get the vertices
				MPI_Win_lock(MPI_LOCK_SHARED, mesh->rankOfVert(j->first), MPI_MODE_NOCHECK, verticesWindow);
				MPI_Get(&localVertices[j->second*3], 3, MPI_DOUBLE,
						mesh->rankOfVert(j->first), mesh->posOfVert(j->first)*3, 3, MPI_DOUBLE,
						verticesWindow);
				MPI_Win_unlock(mesh->rankOfVert(j->first), verticesWindow);
				localVerticesTransfered[j->second] = true;
			}
		}

		// Writing is a collective operation. If we do not have enough partitions,
		// we simply write the last one several times
		unsigned int part = i;
		if (i >= nLocalPart)
			part = nLocalPart-1;

		// Write vertices to the file
		size_t start[3] = {part + rank*nMaxLocalPart, 0, 0};
		checkNcError(nc_put_var1_uint(ncFile, ncVarVrtxSize, start, &localPartVertexSize[part]));
		size_t count[3] = {1, localPartVertexSize[part], 3};
		checkNcError(nc_put_vara_double(ncFile, ncVarVrtxCoords, start, count, localVertices));
	}

	// Remove unused memory
	delete [] localVertices;
	delete [] localVerticesTransfered;

	MPI_Win_free(&verticesWindow);
	MPI_Free_mem(vertices);
	delete [] localPartVerticesMaps;
	delete [] localPartVertexSize;

	logInfo(rank) << "Writing MPI boundary information";
	// Write number of bondaries per partition first
	for (unsigned int i = 0; i < nMaxLocalPart; i++) {
		unsigned int j = i % nLocalPart;

		size_t start[1] = {j + rank*nMaxLocalPart};
		int size = boundaryMaps[j].size();
		checkNcError(nc_put_var1_int(ncFile, ncVarBndSize, start, &size));
	}

	// Count max number of boundaries we have to write on each rank
	maxBoundaries = 0;
	for (unsigned int i = 0; i < nLocalPart; i++) {
		maxBoundaries += boundaryMaps[i].size();
	}
	MPI_Allreduce(MPI_IN_PLACE, &maxBoundaries, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);

	unsigned int* localBoundaryIds = new unsigned int[maxBoundarySize];
	unsigned int boundariesDone = 0;
	for (unsigned int i = 0; i < nLocalPart; i++) {
		unsigned int curBoundary = 0;

		for (std::map<unsigned int, std::vector<MPINeighborElement> >::const_iterator j = boundaryMaps[i].begin();
				j != boundaryMaps[i].end(); j++) {

			size_t start[3] = {i + rank*nMaxLocalPart, curBoundary, 0};
			int size = j->second.size();
			checkNcError(nc_put_var1_int(ncFile, ncVarBndElemSize, start, &size));
			int remoteRank = j->first;
			checkNcError(nc_put_var1_int(ncFile, ncVarBndElemRank, start, &remoteRank));

			// Fill local ids
			for (size_t k = 0; k < j->second.size(); k++) {
				localBoundaryIds[k] = j->second[k].localElement;
			}

			size_t count[3] = {1, 1, static_cast<size_t>(size)};
			checkNcError(nc_put_vara_uint(ncFile, ncVarBndElemLocalIds, start, count, localBoundaryIds));

			boundariesDone++;
			curBoundary++;
		}
	}

	// For collective I/O
	for (size_t k = 0; k < boundaryMaps[0].begin()->second.size(); k++) {
		localBoundaryIds[k] = boundaryMaps[0].begin()->second[k].localElement;
	}

	for (unsigned int i = boundariesDone; i < maxBoundaries; i++) {
		size_t start[3] = {0 + rank*nMaxLocalPart, 0, 0};
		int size = boundaryMaps[0].begin()->second.size();
		checkNcError(nc_put_var1_int(ncFile, ncVarBndElemSize, start, &size));
		int remoteRank = boundaryMaps[0].begin()->first;
		checkNcError(nc_put_var1_int(ncFile, ncVarBndElemRank, start, &remoteRank));

		size_t count[3] = {1, 1, static_cast<size_t>(size)};
		checkNcError(nc_put_vara_uint(ncFile, ncVarBndElemLocalIds, start, count, localBoundaryIds));
	}

	delete [] localBoundaryIds;
	delete [] boundaryMaps;

	delete mesh;

	checkNcError(nc_close(ncFile));

	logInfo(rank) << "Finished successfully";

	PCU_Comm_Free();

	MPI_Finalize();
	return 0;
}
