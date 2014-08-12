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

#ifndef MESH_GENERATOR_H
#define MESH_GENERATOR_H

#include <mpi.h>

#include "MeshInput.h"

/**
 * Superclass for all mesh sources that are based on mesh
 * generators
 */
class MeshGenerator : public MeshInput
{
private:
	/** Number of processes */
	int m_nProcs;

	/** The global id of the first vertex on each process */
	unsigned int* m_vertStart;
	/** The global id of the first element on each process */
	unsigned int* m_elemStart;

public:
	MeshGenerator()
		: m_nProcs(0), m_vertStart(0L), m_elemStart(0L)
	{
	}

	virtual ~MeshGenerator()
	{
		delete [] m_vertStart;
		delete [] m_elemStart;
	}

public:
	unsigned int vertStart(int rank) const
	{
		return m_vertStart[rank];
	}

	int rankOfVert(unsigned int vertex) const
	{
		return std::upper_bound(m_vertStart, m_vertStart+m_nProcs, vertex) - m_vertStart - 1;
	}

	unsigned int posOfVert(unsigned int vertex) const
	{
		return vertex - m_vertStart[rankOfVert(vertex)];
	}

	unsigned int elemStart(int rank) const
	{
		return m_elemStart[rank];
	}

	int rankOfElem(unsigned int element) const
	{
		return std::upper_bound(m_elemStart, m_elemStart+m_nProcs, element) - m_elemStart - 1;
	}

	unsigned int posOfElem(unsigned int element) const
	{
		return element - m_elemStart[rankOfElem(element)];
	}

protected:
	/**
	 * Initializes some additional data structures.
	 * Call this after setting nLocalVertices and nLocalElements
	 */
	void init()
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &m_nProcs);

        // Compute id of the first element for each process
        unsigned int start[2] = {nLocalVertices(), nLocalElements()};
        MPI_Scan(MPI_IN_PLACE, start, 2, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

        m_vertStart = new unsigned int[m_nProcs];
        m_vertStart[rank] = start[0] - nLocalVertices();
        MPI_Allgather(MPI_IN_PLACE, 1, MPI_UNSIGNED, m_vertStart, 1, MPI_UNSIGNED,
        		MPI_COMM_WORLD);

        m_elemStart = new unsigned int[m_nProcs];
        m_elemStart[rank] = start[1] - nLocalElements();
        MPI_Allgather(MPI_IN_PLACE, 1, MPI_UNSIGNED, m_elemStart, 1, MPI_UNSIGNED,
        		MPI_COMM_WORLD);
	}
};

#endif // MESH_GENERATOR_H
