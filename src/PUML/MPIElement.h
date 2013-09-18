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

#ifndef PUML_MPI_ELEMENT_H
#define PUML_MPI_ELEMENT_H

#ifdef PARALLEL
#include <mpi.h>
#endif // PARALLEL

namespace PUML
{

class MPIElement
{
private:
#ifdef PARALLEL
	MPI_Comm m_comm;
#endif // PARALLEL

	int m_rank;

	int m_size;

public:
	MPIElement() :
#ifdef PARALLEL
		m_comm(MPI_COMM_NULL),
#endif // PARALLEL
		m_rank(0), m_size(1)
	{
	}

#ifdef PARALLEL
	void setMPIComm(MPI_Comm comm)
	{
		m_comm = comm;
		MPI_Comm_rank(comm, &m_rank);
		MPI_Comm_size(comm, &m_size);
	}

	MPI_Comm mpiComm()
	{
		return m_comm;
	}
#endif // PARALLEL

	int mpiRank()
	{
		return m_rank;
	}

	int mpiSize()
	{
		return m_size;
	}
};

}

#endif // PUML_MPI_ELEMENT_H
