/**
 * @file
 *  This file is part of PUML
 *
 *  For conditions of distribution and use, please see the copyright
 *  notice in the file 'COPYING' at the root directory of this package
 *  and the copyright notice at https://github.com/TUM-I5/PUML
 *
 * @copyright 2014 Technische Universitaet Muenchen
 * @author Sebastian Rettenberger <rettenbs@in.tum.de>
 */

#ifndef PARALLEL_GAMBIT_READER_H
#define PARALLEL_GAMBIT_READER_H

#include <mpi.h>

#include <vector>

#include "GambitReader.h"

class ParallelGambitReader
{
private:
	const MPI_Comm m_comm;

	int m_rank;
	int m_nProcs;

	GambitReader m_serialReader;

	// Some variables that are required by all processes
	unsigned int m_nVertices;
	unsigned int m_nElements;
	unsigned int m_nBoundaries;

public:
	ParallelGambitReader(MPI_Comm comm = MPI_COMM_WORLD)
		: m_comm(comm), m_nVertices(0), m_nElements(0), m_nBoundaries(0)
	{
		init();
	}

	/**
	 * @param meshFile Only required by rank 0
	 */
	ParallelGambitReader(const char* meshFile, MPI_Comm comm = MPI_COMM_WORLD)
		: m_comm(comm), m_nVertices(0), m_nElements(0), m_nBoundaries(0)
	{
		init();
		open(meshFile);
	}

	void open(const char* meshFile)
	{
		unsigned int vars[3];

		if (m_rank == 0) {
			m_serialReader.open(meshFile);

			vars[0] = m_serialReader.nVertices();
			vars[1] = m_serialReader.nElements();
			vars[2] = m_serialReader.nBoundaries();
		}

		MPI_Bcast(vars, 3, MPI_UNSIGNED, 0, m_comm);

		m_nVertices = vars[0];
		m_nElements = vars[1];
		m_nBoundaries = vars[2];
	}

	unsigned int nVertices() const
	{
		return m_nVertices;
	}

	unsigned int nElements() const
	{
		return m_nElements;
	}

	/**
	 * @return Number of boundary faces
	 */
	unsigned int nBoundaries() const
	{
		return m_nBoundaries;
	}

	/**
	 * Reads all vertices
	 * Each process gets <code>nVertices() + processes - 1) / processes</code> elements,
	 * except for the last, which gets the remaining elements
	 *
	 * This is a collective operation
	 *
	 * @todo Only 3 dimensional meshes are supported
	 */
	void readVertices(double* vertices)
	{
		unsigned int chunkSize = (m_nVertices + m_nProcs - 1) / m_nProcs;

		if (m_rank == 0) {
			// Allocate second buffer so we can read and send in parallel
			double* vertices2 = new double[chunkSize * 3];
			if (m_nProcs % 2 == 0)
				// swap once so we have the correct buffer at the end
				swap(vertices, vertices2);

			MPI_Request request = MPI_REQUEST_NULL;
			for (int i = 1; i < m_nProcs-1; i++) {
				logInfo() << "Reading vertices part" << i << "of" << m_nProcs;
				m_serialReader.readVertices(i * chunkSize, chunkSize, vertices);
				waitIfRequest(request);
				MPI_Isend(vertices, chunkSize*3, MPI_DOUBLE, i, 0, m_comm, &request);
				swap(vertices, vertices);
			}

			// Read last one
			unsigned int lastChunkSize = m_nVertices - (m_nProcs-1) * chunkSize;
			logInfo() << "Reading vertices part" << (m_nProcs-1) << "of" << m_nProcs;
			m_serialReader.readVertices((m_nProcs-1) * chunkSize, lastChunkSize, vertices);
			waitIfRequest(request);
			MPI_Isend(vertices, lastChunkSize*3, MPI_DOUBLE, m_nProcs-1, 0, m_comm, &request);
			swap(vertices, vertices2);

			// Finally read the first part
			logInfo() << "Reading vertices part" << m_nProcs << "of" << m_nProcs;
			m_serialReader.readVertices(0, chunkSize, vertices);
			waitIfRequest(request);

			delete [] vertices2;
		} else {
			if (m_rank == m_nProcs-1)
				chunkSize = m_nVertices - (m_nProcs-1) * chunkSize;

			MPI_Recv(vertices, chunkSize*4, MPI_UNSIGNED, 0, 0, m_comm, MPI_STATUS_IGNORE);
		}
	}


	/**
	 * Reads all elements.
	 * Each process gets <code>nElements() + processes - 1) / processes</code> elements,
	 * except for the last, which gets the remaining elements
	 *
	 * This is a collective operation.
	 *
	 * @todo Only tetrahedral meshes are supported
	 */
	void readElements(unsigned int* elements)
	{
		unsigned int chunkSize = (m_nElements + m_nProcs - 1) / m_nProcs;

		if (m_rank == 0) {
			// Allocate second buffer so we can read and send in parallel
			unsigned int* elements2 = new unsigned int[chunkSize * 4];
			if (m_nProcs % 2 == 0)
				// swap once so we have the correct buffer at the end
				swap(elements, elements2);

			MPI_Request request = MPI_REQUEST_NULL;

			for (int i = 1; i < m_nProcs-1; i++) {
				logInfo() << "Reading elements part" << i << "of" << m_nProcs;
				m_serialReader.readElements(i * chunkSize, chunkSize, elements);
				waitIfRequest(request);
				MPI_Isend(elements, chunkSize*4, MPI_UNSIGNED, i, 0, m_comm, &request);
				swap(elements, elements2);
			}

			// Read last one
			unsigned int lastChunkSize = m_nElements - (m_nProcs-1) * chunkSize;
			logInfo() << "Reading elements part" << (m_nProcs-1) << "of" << m_nProcs;
			m_serialReader.readElements((m_nProcs-1) * chunkSize, lastChunkSize, elements);
			waitIfRequest(request);
			MPI_Isend(elements, lastChunkSize*4, MPI_UNSIGNED, m_nProcs-1, 0, m_comm, &request);
			swap(elements, elements2);

			// Finally read the first part
			logInfo() << "Reading elements part" << m_nProcs << "of" << m_nProcs;
			m_serialReader.readElements(0, chunkSize, elements);
			waitIfRequest(request);

			delete [] elements2;
		} else {
			if (m_rank == m_nProcs-1)
				chunkSize = m_nElements - (m_nProcs-1) * chunkSize;

			MPI_Recv(elements, chunkSize*4, MPI_UNSIGNED, 0, 0, m_comm, MPI_STATUS_IGNORE);
		}
	}


	/**
	 * Reads all group numbers.
	 * Each process gets <code>nElements() + processes - 1) / processes</code> group numbers,
	 * except for the last, which gets the remaining
	 *
	 * This is a collective operation.
	 */
	void readGroups(unsigned int* groups)
	{
		unsigned int chunkSize = (m_nElements + m_nProcs - 1) / m_nProcs;

		if (m_rank == 0) {
			unsigned int maxChunkSize = chunkSize;
			ElementGroup* map = new ElementGroup[maxChunkSize];

			std::vector<unsigned int>* aggregator = new std::vector<unsigned int>[m_nProcs-1];
			unsigned int* sizes = new unsigned int[m_nProcs-1];
			MPI_Request* requests = new MPI_Request[(m_nProcs-1)*2];
			for (int i = 0; i < m_nProcs-1; i++) {
				requests[i*2] = MPI_REQUEST_NULL;
				requests[i*2+1] = MPI_REQUEST_NULL;
			}

			for (int i = 0; i < m_nProcs; i++) {
				logInfo() << "Reading group information part" << (i+1) << "of" << m_nProcs;

				if (i == m_nProcs-1)
					chunkSize = m_nElements - (m_nProcs-1) * chunkSize;

				m_serialReader.readGroups(i*maxChunkSize, chunkSize, map);

				// Wait for all sending from last iteration
				for (int j = 0; j < m_nProcs-1; j++) {
					waitIfRequest(requests[j*2]);
					waitIfRequest(requests[j*2+1]);
					aggregator[j].clear();
				}

				// Sort group numbers into the corresponding aggregator
				for (unsigned int j = 0; j < chunkSize; j++) {
					if (map[j].element < maxChunkSize)
						// Local element
						groups[map[j].element] = map[j].group;
					else {
						// Element for another processor
						// Serials the struct to make sending easier
						unsigned int proc = map[j].element / maxChunkSize;
						aggregator[proc-1].push_back(map[j].element % maxChunkSize);
						aggregator[proc-1].push_back(map[j].group);
					}
				}

				// Send send aggregated mapping
				for (int j = 0; j < m_nProcs-1; j++) {
					if (aggregator[j].empty())
						continue;

					sizes[j] = aggregator[j].size() / 2; // element id + group number
					MPI_Isend(&sizes[j], 1, MPI_UNSIGNED, j+1, 0, m_comm, &requests[j*2]);
					MPI_Isend(&aggregator[j][0], aggregator[j].size(), MPI_UNSIGNED, j+1,
							0, m_comm, &requests[j*2+1]);
				}
			}

			for (int i = 0; i < m_nProcs-1; i++) {
				waitIfRequest(requests[i*2]);
				waitIfRequest(requests[i*2+1]);
			}

			delete [] map;
			delete [] aggregator;
			delete [] sizes;
			delete [] requests;
		} else {
			if (m_rank == m_nProcs-1)
				chunkSize = m_nElements - (m_nProcs-1) * chunkSize;

			// Allocate enough memory
			unsigned int* buf = new unsigned int[chunkSize*2];

			unsigned int recieved = 0;

			while (recieved < chunkSize) {
				unsigned int size;
				MPI_Recv(&size, 1, MPI_UNSIGNED, 0, 0, m_comm, MPI_STATUS_IGNORE);
				MPI_Recv(buf, size*2, MPI_UNSIGNED, 0, 0, m_comm, MPI_STATUS_IGNORE);

				for (unsigned int i = 0; i < size; i += 2)
					groups[buf[i]] = buf[i+1];

				recieved += size;
			}

			delete [] buf;
		}
	}

	/**
	 * Reads all boundaries.
	 * Each process gets the boundaries for <code>nElements() + processes - 1) / processes</code>
	 * elements. Boundaries not specified in the mesh are not modified.
	 *
	 * This is a collective operation
	 *
	 * @todo Only tetrahedral meshes are supported
	 */
	void readBoundaries(unsigned int* boundaries)
	{
		unsigned int chunkSize = (m_nElements + m_nProcs - 1) / m_nProcs;

		if (m_rank == 0) {
			unsigned int maxChunkSize = chunkSize;
			BoundaryFace* faces = new BoundaryFace[maxChunkSize];

			std::vector<unsigned int>* aggregator = new std::vector<unsigned int>[m_nProcs-1];
			unsigned int* sizes = new unsigned int[m_nProcs-1];
			MPI_Request* requests = new MPI_Request[(m_nProcs-1)*2];
			for (int i = 0; i < m_nProcs-1; i++) {
				requests[i*2] = MPI_REQUEST_NULL;
				requests[i*2+1] = MPI_REQUEST_NULL;
			}

			unsigned int nChunks = (m_nBoundaries + chunkSize - 1) / chunkSize;
			for (unsigned int i = 0; i < nChunks; i++) {
				logInfo() << "Reading boundary conditions part" << (i+1) << "of" << nChunks;

				if (i == nChunks-1)
					chunkSize = m_nBoundaries - (nChunks-1) * chunkSize;

				m_serialReader.readBoundaries(i*maxChunkSize, chunkSize, faces);

				// Wait for all sending from last iteration
				for (int j = 0; j < m_nProcs-1; j++) {
					waitIfRequest(requests[j*2]);
					waitIfRequest(requests[j*2+1]);
					aggregator[j].clear();
				}

				// Sort boundary conditions into the corresponding aggregator
				for (unsigned int j = 0; j < chunkSize; j++) {
					if (faces[j].element < maxChunkSize)
						// Local element
						boundaries[faces[j].element*4 + faces[j].face] = faces[j].type;
					else {
						// Face for another processor
						// Serials the struct to make sending easier
						unsigned int proc = faces[j].element / maxChunkSize;
						aggregator[proc-1].push_back((faces[j].element % maxChunkSize)*4
								+ faces[j].face);
						aggregator[proc-1].push_back(faces[j].type);
					}
				}

				// Send send aggregated values
				for (int j = 0; j < m_nProcs-1; j++) {
					if (aggregator[j].empty())
						continue;

					sizes[j] = aggregator[j].size() / 2; // element id + face type
					MPI_Isend(&sizes[j], 1, MPI_UNSIGNED, j+1, 0, m_comm, &requests[j*2]);
					MPI_Isend(&aggregator[j][0], aggregator[j].size(), MPI_UNSIGNED, j+1,
							0, m_comm, &requests[j*2+1]);
				}
			}

			for (int i = 0; i < m_nProcs-1; i++) {
				waitIfRequest(requests[i*2]);
				waitIfRequest(requests[i*2+1]);
			}

			delete [] faces;
			delete [] aggregator;

			// Send finish signal to all other processes
			memset(sizes, 0, (m_nProcs-1)*sizeof(unsigned int));
			for (int i = 0; i < m_nProcs-1; i++)
				MPI_Isend(&sizes[i], 1, MPI_UNSIGNED, i+1, 0, m_comm, &requests[i]);
			MPI_Waitall(m_nProcs-1, requests, MPI_STATUSES_IGNORE);

			delete [] sizes;
			delete [] requests;
		} else {
			if (m_rank == m_nProcs-1)
				chunkSize = m_nElements - (m_nProcs-1) * chunkSize;

			// Allocate enough memory
			unsigned int* buf = new unsigned int[chunkSize*2];

			while (true) {
				unsigned int size;
				MPI_Recv(&size, 1, MPI_UNSIGNED, 0, 0, m_comm, MPI_STATUS_IGNORE);
				if (size == 0)
					// Finished
					break;

				MPI_Recv(buf, size*2, MPI_UNSIGNED, 0, 0, m_comm, MPI_STATUS_IGNORE);

				for (unsigned int i = 0; i < size; i += 2)
					boundaries[buf[i]] = buf[i+1];
			}

			delete [] buf;
		}
	}

private:
	/**
	 * Initialize some parameters
	 */
	void init()
	{
		MPI_Comm_rank(m_comm, &m_rank);
		MPI_Comm_size(m_comm, &m_nProcs);
	}

	/**
	 * Swaps two pointers
	 */
	template<typename T>
	static void swap(T* p1, T* p2)
	{
		T* tmp = p1;
		p1 = p2;
		p2 = tmp;
	}

	/**
	 * Waits for the MPI_Request if it is not NULL
	 */
	static void waitIfRequest(MPI_Request &request)
	{
		if (request != MPI_REQUEST_NULL)
			MPI_Wait(&request, MPI_STATUS_IGNORE);
	}
};

#endif // PARALLEL_GAMBIT_READER_H
