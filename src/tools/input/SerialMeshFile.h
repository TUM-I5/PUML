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

#ifndef SERIAL_MESH_FILE_H
#define SERIAL_MESH_FILE_H

#ifdef PARALLEL
#include <mpi.h>
#endif // PARALLEL

/**
 * Read a mesh from a serial file
 */
template<typename T>
class SerialMeshFile
{
private:
#ifdef PARALLEL
	MPI_Comm m_comm;
#endif // PARALLEL

	int m_rank;
	int m_nProcs;

	T m_meshReader;

	/** Number of local vertices */
	unsigned int m_nLocalVertices;
	/** Maximum number of local vertices */
	unsigned int m_nMaxLocalVertices;
	/** Number of local elements */
	unsigned int m_nLocalElements;
	/** Maximum number of local elements */
	unsigned int m_nMaxLocalElements;

public:
#ifdef PARALLEL
	SerialMeshFile(MPI_Comm comm = MPI_COMM_WORLD)
		: m_comm(comm), m_meshReader(comm)
	{
		init();
	}
#else // PARALLEL
	SerialMeshFile()
		: m_nLocalElements(0)
	{
		init();
	}
#endif // PARALLEL

#ifdef PARALLEL
	SerialMeshFile(const char* meshFile, MPI_Comm comm = MPI_COMM_WORLD)
		: m_comm(comm), m_meshReader(comm)
	{
		init();
		open(meshFile);
	}
#else // PARALLEL
	SerialMeshFile(const char* meshFile)
	{
		init();
		open(meshFile);
	}
#endif // PARALLEL

	void open(const char* meshFile)
	{
		m_meshReader.open(meshFile);

		m_nMaxLocalVertices = (nVertices + m_nProcs - 1) / m_nProcs;
		m_nMaxLocalElements = (nElements() + m_nProcs - 1) / m_nProcs;
		if (m_rank == m_nProcs - 1) {
			m_nLocalVertices = nVertices() - (m_nProcs-1) * m_nMaxLocalVertices;
			m_nLocalElements = nElements() - (m_nProcs-1) * m_nMaxLocalElements;
		} else {
			m_nLocalVertices = m_nMaxLocalVertices;
			m_nLocalElements = m_nMaxLocalElements;
		}
	}

	unsigned int nVertices() const
	{
		return m_meshReader.nVertices();
	}

	unsigned int nElements() const
	{
		return m_meshReader.nElements();
	}

	/**
	 * Number of vertices this process is responsible for
	 */
	unsigned int nLocalVertices() const
	{
		return m_nLocalVertices;
	}

	/**
	 * Number of elements this process is responsible for
	 */
	unsigned int nLocalElements() const
	{
		return m_nLocalElements;
	}

	/**
	 * @return The rank of this vertex
	 */
	int rankOfVert(unsigned int vertex) const
	{
		return vertex / m_nMaxLocalElements;
	}

	/**
	 * @return The position of the vertex on the local rank
	 */
	unsigned int posOfVert(unsigned int vertex) const
	{
		return vertex % m_nMaxLocalElements;
	}

	/**
	 * @return The first element id for process <code>rank</code>
	 */
	unsigned int elemStart(int rank) const
	{
		return rank * m_nMaxLocalElements;
	}

	/**
	 * @return The rank for this element
	 */
	int rankOfElem(unsigned int element) const
	{
		return element / m_nMaxLocalElements;
	}

	/**
	 * @return The position of the element on the local rank
	 */
	unsigned int posOfElem(unsigned int element) const
	{
		return element % m_nMaxLocalElements;
	}

	/**
	 * Get vertices
	 */
	void getVertices(double* vertices)
	{
		m_meshReader.readVertices(vertices);
	}

	/**
	 * Get elements
	 */
	void getElements(unsigned int* elements)
	{
		m_meshReader.readElements(elements);
	}

	/**
	 * Get group information for the elements
	 */
	void getGroups(unsigned int* groups)
	{
		m_meshReader.readGroups(groups);
	}

	/**
	 * Get boundary information
	 */
	void getBoundaries(unsigned int* boundaries)
	{
		m_meshReader.readBoundaries(boundaries);
	}

private:
	/**
	 * Sets some parameters (called from the constructor)
	 */
	void init()
	{
#ifdef PARALLEL
		MPI_Comm_rank(m_comm, &m_rank);
		MPI_Comm_size(m_comm, &m_nProcs);
#else // PARALLLEL
		m_rank = 0;
		m_nProcs = 1;
#endif // PARALLEL
	}
};

#endif // SERIAL_MESH_FILE_H
