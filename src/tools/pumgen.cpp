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

#ifdef PARALLEL
#include <mpi.h>
#endif // PARALLEL

#include <algorithm>
#include <cstring>
#include <map>
#include <string>
#include <vector>

#include "utils/args.h"
#include "utils/logger.h"

#ifdef PARALLEL
#include <parmetis.h>
#else // PARALLEL
#include <metis.h>
#endif // PARALLEL

#ifdef PARALLEL
#include <netcdf_par.h>
#endif // PARALLEL
#include <netcdf.h>

#include "meshreader/GambitReader.h"

struct Element {
	/** Global id */
	unsigned int id;
	/** Vertex ids */
	unsigned int vertex[4];
	/** The group number */
	unsigned int group;

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

int main(int argc, char* argv[])
{
	int rank = 0;
	int processes = 1;

#ifdef PARALLEL
	int threadSupport;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &threadSupport);

	if (threadSupport != MPI_THREAD_SINGLE)
		// I'm not sure yet if we really need thread support
		logError() << "Threading is not supported in this MPI version";

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &processes);
#endif // PARALLEL

	GambitReader gambitReader;
	unsigned int nPartitions;
	std::string outputFile;

	unsigned int nVertices;
	unsigned int nElements;

	unsigned int nLocalElements;	// The number of elements on this rank
	unsigned int nMaxLocalElements;	// The number of elements on all ranks except the last

	// TODO only tetrahedral meshes are currently supported
	unsigned int* elements;

	if (rank == 0) {
		std::string inputFile;

		utils::Args args;
		args.addAdditionalOption("input", "Input mesh file");
		args.addAdditionalOption("partition", "Number of partitions");
		args.addAdditionalOption("output", "Output parallel unstructured mesh file", false);

		if (args.parse(argc, argv) == utils::Args::Success) {
			inputFile = args.getAdditionalArgument<std::string>("input");
			nPartitions = args.getAdditionalArgument<unsigned int>("partition");

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

			if (nPartitions == 0)
				logError() << "Partitions created must be greater than zero";

			if (nPartitions <= ((nPartitions + processes - 1) / processes) * (processes-1))
				logError() << "Not every process will get at least one partition, use a smaller number of ranks";

#ifdef PARALLEL
			int constants[] = {nPartitions, static_cast<int>(outputFile.size())};
			MPI_Bcast(constants, 2, MPI_INT, 0, MPI_COMM_WORLD);

			char* buf = new char[outputFile.size()];
			outputFile.copy(buf, outputFile.size());
			MPI_Bcast(buf, outputFile.size(), MPI_CHAR, 0, MPI_COMM_WORLD);
			delete [] buf;
#endif // PARALLEL
		} else {
#ifdef PARALLEL
			// Tell all processes to quit
			int constants[] = {0, 0};
			MPI_Bcast(constants, 1, MPI_INT, 0, MPI_COMM_WORLD);

			MPI_Finalize();
#endif // PARALLEL
			return 1;
		}

		// Rank 0 does all the reading and distributes the contents
		gambitReader.open(inputFile.c_str());

		nVertices = gambitReader.nVertices();
		nElements = gambitReader.nElements();

#ifdef PARALLEL
		// Distribute the mesh size
		unsigned int meshSize[] = {nVertices, nElements};
		MPI_Bcast(meshSize, 2, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
#endif // PARALLEL

		nMaxLocalElements = (nElements + processes - 1) / processes;
		nLocalElements = nMaxLocalElements;

		// Allocate memory for the elements
#ifdef PARALLEL
		MPI_Alloc_mem(sizeof(unsigned int)*nLocalElements*4, MPI_INFO_NULL, &elements);
#else // PARALLEL
		elements = new unsigned int[nLocalElements * 4];
#endif // PARALLEL

#ifdef PARALLEL
		// Read chunks of elements and distribute them to each processor
		for (int i = 1; i < processes-1; i++) {
			logInfo() << "Reading elements part" << i << "of" << processes;
			gambitReader.readElements(i * nLocalElements, nLocalElements, elements);
			MPI_Send(elements, nLocalElements*4, MPI_UNSIGNED, i, 0, MPI_COMM_WORLD);
		}

		// Last process
		if (processes > 1) {
			logInfo() << "Reading elements part" << (processes-1) << "of" << processes;
			gambitReader.readElements((processes-1) * nLocalElements, nElements - (processes-1) * nLocalElements, elements);
			MPI_Send(elements, (nElements - (processes-1) * nLocalElements) * 4, MPI_UNSIGNED, processes - 1, 0, MPI_COMM_WORLD);
		}
#endif // PARALLEL

		// Finally our own elements
		logInfo() << "Reading elements part" << processes << "of" << processes;
		gambitReader.readElements(0, nLocalElements, elements);
	} else {
#ifdef PARALLEL
		// All ranks except 0
		int constants[2];
		MPI_Bcast(constants, 2, MPI_INT, 0, MPI_COMM_WORLD);

		if (constants[0] == 0) {
			// Wrong command line parameters
			MPI_Finalize();
			return 1;
		}

		nPartitions = constants[0];

		char* buf = new char[constants[1]];
		MPI_Bcast(buf, constants[1], MPI_CHAR, 0, MPI_COMM_WORLD);
		outputFile.append(buf, constants[1]);
		delete [] buf;

		// Size of the mesh
		unsigned int meshSize[2];
		MPI_Bcast(meshSize, 2, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
		nVertices = meshSize[0];
		nElements = meshSize[1];

		nMaxLocalElements = (nElements + processes - 1) / processes;
		if (rank == processes - 1)
			nLocalElements = nElements - (processes-1) * nMaxLocalElements;
		else
			nLocalElements = nMaxLocalElements;

		MPI_Alloc_mem(sizeof(unsigned int)*nLocalElements*4, MPI_INFO_NULL, &elements);

		// Get our part from rank 0
		MPI_Recv(elements, nLocalElements*4, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#endif // PARALLEL
	}

	// TODO wrap this code part so we can switch to different libraries

	logInfo(rank) << "Creating dual graph with METIS";
	idx_t* elemDist = new idx_t[processes+1];
	for (int i = 0; i < processes; i++)
		elemDist[i] = i * nMaxLocalElements;
	elemDist[processes] = nElements;

	idx_t* eptr = new idx_t[nLocalElements+1];
	for (unsigned int i = 0; i < nLocalElements+1; i++)
		eptr[i] = i * 4;

	idx_t* eind = new idx_t[nLocalElements*4];
	for (unsigned int i = 0; i < nLocalElements*4; i++)
		eind[i] = elements[i];

	idx_t numflag = 0;
	idx_t ncommonnodes = 3;

	idx_t* xadj;
	idx_t* adjncy;
#ifdef PARALLEL
	MPI_Comm commWorld = MPI_COMM_WORLD;
	if (ParMETIS_V3_Mesh2Dual(elemDist, eptr, eind, &numflag, &ncommonnodes, &xadj, &adjncy, &commWorld) != METIS_OK)
#else // PARALLEL
    idx_t ne = nElements;
	idx_t nn = nVertices;
    if (METIS_MeshToDual(&ne, &nn, eptr, eind, &ncommonnodes, &numflag, &xadj, &adjncy) != METIS_OK)
#endif // PARALLEL
		logError() << "Could not create dual graph";

	delete [] eptr;
	delete [] eind;

	logInfo(rank) << "Creating neighborhood information";
	int* elementNeighbors;
#ifdef PARALLEL
	MPI_Alloc_mem(sizeof(unsigned int)*4*nLocalElements, MPI_INFO_NULL, &elementNeighbors);
#else // PARALLEL
	elementNeighbors = new unsigned int[4*nLocalElements];
#endif // PARALLEL))
	for (unsigned int i = 0; i < nLocalElements; i++) {
		idx_t j;
		for (j = 0; j < xadj[i+1]-xadj[i]; j++)
			elementNeighbors[i*4 + j] = adjncy[xadj[i] + j];

		for (; j < 4; j++)
			elementNeighbors[i*4 + j] = -1;

	}
#ifdef PARALLEL
	MPI_Win elementNbWindow;
	MPI_Win_create(elementNeighbors, sizeof(unsigned int)*nLocalElements*4,
			sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &elementNbWindow);
#endif // PARALLEL

	logInfo(rank) << "Creating partitions with METIS";
	idx_t wgtflag = 0;
	idx_t ncon = 1;
	idx_t nparts = nPartitions;

	real_t* tpwgts = new real_t[nPartitions];
	for (unsigned int i = 0; i < nPartitions; i++)
		tpwgts[i] = 1./nPartitions;

	real_t ubev = 1.01;
	idx_t options = 0;

	idx_t edgecut;
	idx_t* part = new idx_t[nLocalElements];

#ifdef PARALLEL
	if (ParMETIS_V3_PartKway(elemDist, xadj, adjncy, 0L, 0L, &wgtflag, &numflag, &ncon,
			&nparts, tpwgts, &ubev, &options, &edgecut, part, &commWorld) != METIS_OK)
#else // PARALLEL
	// TODO
#endif // PARALLEL
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

	// Read and distribute group information
	unsigned int* elementGroups;
#ifdef PARALLEL
	MPI_Alloc_mem(sizeof(unsigned int)*nLocalElements, MPI_INFO_NULL, &elementGroups);

	MPI_Win groupWindow;
	MPI_Win_create(elementGroups, sizeof(unsigned int)*nLocalElements, sizeof(unsigned int),
			MPI_INFO_NULL, MPI_COMM_WORLD, &groupWindow);
#else // PARALLEL
	elementGroups = new unsigned int[nLocalElements];
#endif // PARALLEL


	if (rank == 0) {
#ifdef PARALLEL
		ElementGroup* groups = new ElementGroup[nLocalElements];

		// Read the group and put the on the other processes
		for (int i = 0; i < processes; i++) {
			unsigned int count = nLocalElements;
			if (i == processes-1)
				count = nElements - (processes-1) * nLocalElements;

			logInfo() << "Reading group information part" << (i+1) << "of" << processes;
			gambitReader.readGroups(i * nLocalElements, count, groups);
			for (unsigned int j = 0; j < count; j++) {
				MPI_Win_lock(MPI_LOCK_EXCLUSIVE, groups[j].element / nLocalElements, MPI_MODE_NOCHECK, groupWindow);
				MPI_Put(&groups[j].group, 1, MPI_UNSIGNED, groups[j].element / nLocalElements, groups[j].element % nLocalElements,
						1, MPI_UNSIGNED, groupWindow);
				MPI_Win_unlock(groups[j].element / nLocalElements, groupWindow);
			}
		}

		delete [] groups;
#else // PARALLEL
		// TODO
#endif // PARALLEL
	}

#ifdef PARALLEL
	MPI_Win_free(&groupWindow);
#endif // PARALLEL

	logInfo(rank) << "Redistributing elements";
	unsigned int nMaxLocalPart = (nPartitions + processes - 1) / processes;
	unsigned int nLocalPart = nMaxLocalPart;
#ifdef PARALLEL
	if (rank == processes - 1)
		nLocalPart = nPartitions - (processes-1) * nMaxLocalPart;
#endif // PARALLEL

	// Redistribute the elements among the processors according to there partitions
	// Count number local elements in each partition
	unsigned int* localPartPtr = new unsigned int[nLocalPart+1];
	localPartPtr[0] = 0;

	Element* localPartElements;

	unsigned int pos;	// Helper variable

#ifdef PARALLEL
	std::vector<unsigned int>* elementPerPart = new std::vector<unsigned int>[nPartitions];
	for (unsigned int i = 0; i < nLocalElements; i++)
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

	// Now send the ids+vertices of the elements to all processes
	Element* sendElements = new Element[nLocalElements];
	pos = 0;
	for (unsigned int i = 0; i < nPartitions; i++) {
		for (unsigned int j = 0; j < nElementPerPart[i]; j++) {
			// Set global id
			sendElements[pos].id = elementPerPart[i][j] + nMaxLocalElements * rank;
			memcpy(sendElements[pos].vertex, &elements[elementPerPart[i][j]*4], sizeof(unsigned int)*4);
			sendElements[pos].group = elementGroups[elementPerPart[i][j]];
			pos++;
		}
	}

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
	int blockLength[] = {1, 4, 1, 1};
	MPI_Aint displacement[] = {0, reinterpret_cast<uintptr_t>(recvElements[0].vertex)
			- reinterpret_cast<uintptr_t>(&recvElements[0]),
			reinterpret_cast<uintptr_t>(&recvElements[0].group)
			- reinterpret_cast<uintptr_t>(&recvElements[0]), sizeof(Element)};
	MPI_Datatype type[] = {MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED, MPI_UB};
	MPI_Datatype elementType;
	MPI_Type_create_struct(4, blockLength, displacement, type, &elementType);
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
#ifdef PARALLEL
	MPI_Free_mem(elementGroups);
#else // PARALLEL
	delete [] elementGroups;
#endif // PARALLEL

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
#else // PARALLEL
	for (unsigned int i = 0; i < nLocalElements; i++)
		localPartSize[part[i]]++;
	localSize = nLocalElements;

	localPartElements = new Element[localSize];
	// TODO fill localPartElements
#endif // PARALLEL

	// Save the partition information because other processes still need it
	// But copy it to an unsigned int
	unsigned int* elementPart;
#ifdef PARALLEL
	MPI_Alloc_mem(sizeof(unsigned int)*nLocalElements, MPI_INFO_NULL, &elementPart);
#else // PARALLEL
	elementPart = new unsigned int[nLocalElements];
#endif // PARALLEL
	for (unsigned int i = 0; i < nLocalElements; i++)
		elementPart[i] = part[i];

	delete [] part;

	logInfo(rank) << "Sorting element ids";
	for (unsigned int i = 0; i < nLocalPart; i++)
		std::sort(&localPartElements[localPartPtr[i]],
				&localPartElements[localPartPtr[i]]+localPartPtr[i+1]-localPartPtr[i]);

	// Read boundary information from the mesh
	unsigned int* boundaries;
#ifdef PARALLEL
	MPI_Alloc_mem(sizeof(unsigned int)*nLocalElements*4, MPI_INFO_NULL, &boundaries);
	MPI_Win boundaryWindow;
	MPI_Win_create(boundaries, sizeof(unsigned int)*nLocalElements*4,
			sizeof(unsigned int), MPI_INFO_NULL, MPI_COMM_WORLD, &boundaryWindow);
#else // PARALLEL
	boundaries = new unsigned int[nLocalElements*4];
#endif // PARELLEL
	memset(boundaries, 0, sizeof(unsigned int)*nLocalElements*4);
	if (rank == 0) {
		// Rank 0 does all the reading and distributes the contents

		unsigned int nBoundaries = gambitReader.nBoundaries();

		// Load the boundaries in chunks
		BoundaryFace* boundaryFaces = new BoundaryFace[nLocalElements];

		for (unsigned int i = 0; i < (nBoundaries + nLocalElements - 1) / nLocalElements; i++) {
			logInfo() << "Reading boundary part" << (i+1)
					<< "of" << ((nBoundaries + nLocalElements - 1) / nLocalElements);
			unsigned int n = std::min(nLocalElements, nBoundaries - i*nLocalElements);
			gambitReader.readBoundaries(i * nLocalElements, n, boundaryFaces);

			// Store the boundary information in the correct rank
			for (unsigned int j = 0; j < n; j++) {
				boundaryFaces[j].type -= 100; // TODO do not do this here

#ifdef PARALLEL
				MPI_Win_lock(MPI_LOCK_EXCLUSIVE, boundaryFaces[j].element / nMaxLocalElements,
						MPI_MODE_NOCHECK, boundaryWindow);
				MPI_Put(&boundaryFaces[j].type, 1, MPI_UNSIGNED, boundaryFaces[j].element / nMaxLocalElements,
						(boundaryFaces[j].element % nMaxLocalElements)*4 + boundaryFaces[j].face,
						1, MPI_UNSIGNED, boundaryWindow);
				MPI_Win_unlock(boundaryFaces[j].element / nMaxLocalElements, boundaryWindow);
#else // PARALLEL
				boundaries[(boundaryFaces[j].element % nMaxLocalElements)*4 + boundaryFaces[j].face] = boundaryFaces[j].type;
#endif // PARALLEL
			}
		}

		delete [] boundaryFaces;
	}

	logInfo(rank) << "Computing element data";
#ifdef PARALLEL
	// Create MPI window to allow access to all elements
	MPI_Win elementWindow;
	MPI_Win_create(elements, sizeof(unsigned int)*nLocalElements*4,
			sizeof(unsigned int), MPI_INFO_NULL, MPI_COMM_WORLD, &elementWindow);

	// Create window for partition information
	MPI_Win elementPartWindow;
	MPI_Win_create(elementPart, sizeof(unsigned int)*nLocalElements,
			sizeof(unsigned int), MPI_INFO_NULL, MPI_COMM_WORLD, &elementPartWindow);
#endif // PARALLEL

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

		int neighbors[4];
#ifdef PARALLEL
		MPI_Win_lock(MPI_LOCK_SHARED, localPartElements[i].id / nMaxLocalElements, MPI_MODE_NOCHECK, elementNbWindow);
		MPI_Get(neighbors, 4, MPI_INT,
				localPartElements[i].id / nMaxLocalElements, (localPartElements[i].id % nMaxLocalElements)*4, 4, MPI_INT,
				elementNbWindow);
		MPI_Win_unlock(localPartElements[i].id / nMaxLocalElements, elementNbWindow);

		// Set boundary (do this for all boundaries, because of inner boundary conditions)
		MPI_Win_lock(MPI_LOCK_SHARED, localPartElements[i].id / nMaxLocalElements, MPI_MODE_NOCHECK, boundaryWindow);
		MPI_Get(&localPartElemBoundary[i*4], 4, MPI_UNSIGNED,
				localPartElements[i].id / nMaxLocalElements, (localPartElements[i].id % nMaxLocalElements)*4, 4, MPI_UNSIGNED,
				boundaryWindow);
		MPI_Win_unlock(localPartElements[i].id / nMaxLocalElements, boundaryWindow);
#else // PARALLEL
		// TODO
		memcpy(localPartElemBoundary[i*4], boundaries[localPartIds[i]*4 ], sizeof(unsigned int)*4);
#endif // PARALLEL

		for (int j = 0; j < 4; j++) {
			if (neighbors[j] < 0)
				continue;

			// Check whether we are in the same partition
			unsigned int neighborPart;
#ifdef PARALLEL
			MPI_Win_lock(MPI_LOCK_SHARED, neighbors[j] / nMaxLocalElements, MPI_MODE_NOCHECK, elementPartWindow);
			MPI_Get(&neighborPart, 1, MPI_UNSIGNED,
					neighbors[j] / nMaxLocalElements, neighbors[j] % nMaxLocalElements, 1, MPI_UNSIGNED,
					elementPartWindow);
			MPI_Win_unlock(neighbors[j] / nMaxLocalElements, elementPartWindow);
#else // PARALLEL
			// TODO
#endif // PARALLEL

			unsigned int neighborVertices[4];
#ifdef PARALLEL
			MPI_Win_lock(MPI_LOCK_SHARED, neighbors[j] / nMaxLocalElements, MPI_MODE_NOCHECK, elementWindow);
			MPI_Get(neighborVertices, 4, MPI_UNSIGNED,
					neighbors[j] / nMaxLocalElements, (neighbors[j] % nMaxLocalElements)*4, 4, MPI_UNSIGNED,
					elementWindow);
			MPI_Win_unlock(neighbors[j] / nMaxLocalElements, elementWindow);
#else // PARALLEL
			// TODO
#endif // PARALLEL

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
				localPartElemNb[i*4 + face] = globalElem2PartId.at(neighbors[j]);
			else {
				// Mark id as neighbor
				MPINeighborElement neighborElement = {i-localPartPtr[curPart], face, neighbors[j], face2};
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

#ifdef PARALLEL
	MPI_Win_free(&boundaryWindow);
	MPI_Free_mem(boundaries);
#else // PARALLEL
	delete [] boundaries;
#endif // PARALLEL

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
#ifdef PARALLEL
	MPI_Allreduce(MPI_IN_PLACE, &maxPartSize, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
#endif // PARALLEL

	// Compute max number of vertices in a partition
	unsigned int maxVertices = localPartVertexSize[0];
	for (unsigned int i = 1; i < nLocalPart; i++)
		if (maxVertices < localPartVertexSize[i])
			maxVertices = localPartVertexSize[i];
#ifdef PARALLEL
	MPI_Allreduce(MPI_IN_PLACE, &maxVertices, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
#endif // PARALLEL

	// Compute max number of boundaries per partition
	unsigned int maxBoundaries = boundaryMaps[0].size();
	for (unsigned int i = 1; i < nLocalPart; i++)
		if (maxBoundaries < boundaryMaps[i].size())
			maxBoundaries = boundaryMaps[i].size();
#ifdef PARALLEL
	MPI_Allreduce(MPI_IN_PLACE, &maxBoundaries, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
#endif // PARALLEL

	// Compute max boundary size
	unsigned int maxBoundarySize = 0;
	for (unsigned int i = 0; i < nLocalPart; i++) {
		for (std::map<unsigned int, std::vector<MPINeighborElement> >::const_iterator j = boundaryMaps[i].begin();
				j != boundaryMaps[i].end(); j++)
			if (maxBoundarySize < j->second.size())
				maxBoundarySize = j->second.size();

	}
#ifdef PARALLEL
	MPI_Allreduce(MPI_IN_PLACE, &maxBoundarySize, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
#endif // PARALLEL

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
#ifdef PARALLEL
	checkNcError(nc_create_par(outputFile.c_str(), NC_NETCDF4 | NC_MPIIO, MPI_COMM_WORLD, MPI_INFO_NULL, &ncFile));
#else // PARALLEL
	checkNcError(nc_create(outputFile.c_str(), NC_NETCDF4, &ncFile));
#endif

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
		size_t count[3] = {1, size, 4};
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

#ifdef PARALLEL
	MPI_Win_free(&elementWindow);
	MPI_Win_free(&elementPartWindow);
	MPI_Win_free(&elementNbWindow);
	MPI_Free_mem(elements);
	MPI_Free_mem(elementPart);
	MPI_Free_mem(elementNeighbors);
#else // PARALLEL
	delete [] elements;
	delete [] elementPart;
	delete [] elementNeighbors;
#endif // PARALLEL
	delete [] localPartPtr;
	delete [] localPartElements;

	// Load vertices and distribute them to the other processes
	unsigned int nMaxLocalVertices = (nVertices + processes - 1) / processes;
	unsigned int nLocalVertices = nMaxLocalVertices;
#ifdef PARALLEL
	if (rank == processes - 1)
		nLocalVertices = nVertices - nMaxLocalVertices * (processes - 1);
#endif // PARALLEL

	double* vertices;
#ifdef PARALLEL
		MPI_Alloc_mem(sizeof(double)*nLocalVertices*3, MPI_INFO_NULL, &vertices);
#else // PARALLEL
		vertices = new double[nLocalVertices * 3];
#endif // PARALLEL


	if (rank == 0) {
#ifdef PARALLEL
		// Read the vertices and send them to the other processes
		for (int i = 1; i < processes-1; i++) {
			logInfo() << "Reading vertices part" << i << "of" << processes;
			gambitReader.readVertices(i * nLocalVertices, nLocalVertices, vertices);
			MPI_Send(vertices, nLocalVertices*3, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
		}

		// Last process
		if (processes > 1) {
			logInfo() << "Reading vertices part" << (processes-1) << "of" << processes;
			gambitReader.readVertices((processes-1) * nLocalVertices, nVertices - (processes-1) * nLocalVertices, vertices);
			MPI_Send(vertices, (nVertices - (processes-1) * nLocalVertices) * 3, MPI_DOUBLE, processes - 1, 0, MPI_COMM_WORLD);
		}
#endif // PARALLEL

		// Finally our own vertices
		logInfo() << "Reading vertices part" << processes << "of" << processes;
		gambitReader.readVertices(0, nLocalVertices, vertices);
	} else {
#ifdef PARALLEL
		MPI_Recv(vertices, nLocalVertices*3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#endif // PARALLEL
	}

#ifdef PARALLEL
	// Create MPI window to get access to all vertices
	MPI_Win verticesWindow;
	MPI_Win_create(vertices, sizeof(double)*nLocalVertices,
			sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &verticesWindow);
#endif // PARALLEL

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
#ifdef PARALLEL
				MPI_Win_lock(MPI_LOCK_SHARED, j->first / nMaxLocalVertices, MPI_MODE_NOCHECK, verticesWindow);
				MPI_Get(&localVertices[j->second*3], 3, MPI_DOUBLE,
						j->first / nMaxLocalVertices, (j->first % nMaxLocalVertices)*3, 3, MPI_DOUBLE,
						verticesWindow);
				MPI_Win_unlock(j->first / nMaxLocalVertices, verticesWindow);
#else // PARALLEL
				// TODO
#endif // PARALLEL
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

#ifdef PARALLEL
	MPI_Win_free(&verticesWindow);
	MPI_Free_mem(vertices);
#else // PARALLEL
	delete [] vertices;
#endif // PARALLEL
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
#ifdef PARALLEL
	MPI_Allreduce(MPI_IN_PLACE, &maxBoundaries, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
#endif // PARALLEL

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

			size_t count[3] = {1, 1, size};
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

		size_t count[3] = {1, 1, size};
		checkNcError(nc_put_vara_uint(ncFile, ncVarBndElemLocalIds, start, count, localBoundaryIds));
	}

	delete [] localBoundaryIds;
	delete [] boundaryMaps;

	checkNcError(nc_close(ncFile));

	logInfo(rank) << "Finished succesfully";

#ifdef PARALLEL
	MPI_Finalize();
#endif // PARALLEL
	return 0;
}
