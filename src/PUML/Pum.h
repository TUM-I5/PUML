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

#ifndef PUML_PUM_H
#define PUML_PUM_H

#ifdef PARALLEL
#include <mpi.h>
#endif // PARALLEL

#include <cstdlib>
#include <map>
#include <string>
#include <vector>

#include "PUML/Group.h"
#include "PUML/MPIElement.h"

namespace PUML
{

class Pum : protected MPIElement
{
private:
	/** Number of partitions */
	size_t m_numPartitions;

public:
	Pum()
		: m_numPartitions(0)
	{
	}

	virtual ~Pum()
	{ }

	/**
	 * Create a new PUM file in serial mode
	 *
	 * @param numPartitions The number of partitions in the mesh
	 * @return True on success, false otherwise
	 */
	bool create(const char* path, size_t numPartitions)
	{
		m_numPartitions = numPartitions;
		return _create(path);
	}

#ifdef PARALLEL
	/**
	 * Create a new PUM file in parallel mode
	 *
	 * @param numPartitions The number of partitions in the mesh (may differ from the number of ranks)
	 * @return True on success, false otherwise
	 */
	bool create(const char* path, size_t numPartitions, MPI_Comm comm, MPI_Info info = MPI_INFO_NULL)
	{
		m_numPartitions = numPartitions;
		setMPIComm(comm);
		return _create(path, comm, info);
	}
#endif // PARALLEL

	virtual bool open(const char* path) = 0;

#ifdef PARALLEL
	bool open(const char* path, MPI_Comm comm, MPI_Info info = MPI_INFO_NULL)
	{
		setMPIComm(comm);
		return _open(path, comm, info);
	}
#endif // PARALLEL

	/**
	 * Create a group in this file
	 *
	 * @param size The total size (number of all elements in all partitions) of this group.
	 *  Use Group::UNLIMITED if unknown
	 */
	virtual Group* createGroup(const char* name, size_t size = Group::UNLIMITED) = 0;

	/**
	 * Create an indexed group in this file
	 *
	 * @see createGroup
	 */
	virtual Group* createGroupIndexed(const char* name, size_t size = Group::UNLIMITED, size_t indexSize = Group::UNLIMITED)
	{
		Group* group = createGroup(name, size);
		if (!group)
			return 0L;

		group->addIndex(indexSize);

		return group;
	}

	/**
	 * Creates a special group for vertices
	 *
	 * @ingroup HighLevelApi
	 */
	Group* createVertexGroup(size_t size = Group::UNLIMITED, size_t indexSize = Group::UNLIMITED)
	{
		return createGroupIndexed("vertex", size, indexSize);
	}

	/**
	 * Creates a special group for cells
	 *
	 * @ingroup HighLevelApi
	 */
	Group* createCellGroup(size_t size = Group::UNLIMITED)
	{
		return createGroup("cell", size);
	}

	virtual Group* getGroup(const char* name) = 0;

	/**
	 * End the definition phase. Groups and entities can only be added during the
	 * definition phase
	 */
	virtual bool endDefinition()
	{
		return true;
	}

	virtual bool close() = 0;

	/**
	 * @return Number of partitions in this file
	 */
	size_t numPartitions() const
	{
		return m_numPartitions;
	}

protected:
	void setNumPartitions(size_t numPartitions)
	{
		m_numPartitions = numPartitions;
	}

	virtual bool _create(const char* path) = 0;
#ifdef PARALLEL
	virtual bool _create(const char* path, MPI_Comm comm, MPI_Info info = MPI_INFO_NULL) = 0;
	virtual bool _open(const char* path, MPI_Comm comm, MPI_Info info = MPI_INFO_NULL) = 0;
#endif // PARALLEL

protected:
	static const std::string CONVENTIONS;
	static const int FILE_VERSION;

	static const char* ATT_CONVENTIONS;
	static const char* ATT_FILE_VERSION;
	static const char* ATT_NUM_PARTITIONS;
};

}

#endif // PUML_PUM_H
