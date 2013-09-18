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

#ifndef PUML_NETCDF_GROUP_H
#define PUML_NETCDF_GROUP_H

#include <netcdf.h>

#include "PUML/Dimension.h"
#include "PUML/Group.h"
#include "PUML/NetcdfElement.h"
#include "PUML/NetcdfEntity.h"

namespace PUML
{

class NetcdfPum;

class NetcdfGroup : public Group, public NetcdfElement
{
private:
	/** nc identifier for the partitions dimension */
	int m_ncDimPartition;

	/** nc identifier for the size dimension */
	int m_ncDimSize;

	/** User defined dimensions */
	std::vector<Dimension> m_dimensions;

	/** nc identifier for the offset variable */
	int m_ncVarOffset;

	/** Entities in this group */
	std::vector<NetcdfEntity> m_entities;

public:
	NetcdfGroup()
		: m_ncDimPartition(-1), m_ncDimSize(-1), m_ncVarOffset(-1)
	{
	}

	/**
	 * @param size The total size of this group. Use NC_UNLIMITED if unknown
	 */
	NetcdfGroup(const char* name, size_t numPartitions, NetcdfElement &ncPum, MPIElement comm, size_t size = NC_UNLIMITED)
		: Group(name, numPartitions, comm), NetcdfElement(&ncPum)
	{
		int ncGroup;
		if (checkError(nc_def_grp(ncPum.identifier(), name, &ncGroup)))
			return;

		setIdentifier(ncGroup);

		if (checkError(nc_def_dim(identifier(), DIM_PARTITION, numPartitions, &m_ncDimPartition)))
			return;

		if (checkError(nc_def_dim(identifier(), DIM_SIZE, size, &m_ncDimSize)))
			return;

		if (checkError(nc_def_var(identifier(), VAR_OFFSET, NC_UINT64, 1, &m_ncDimPartition, &m_ncVarOffset)))
			return;
#ifdef PARALLEL
		if (checkError(nc_var_par_access(identifier(), m_ncVarOffset, NC_COLLECTIVE)))
			return;
#endif // PARALLEL
	}

	Dimension& createDimension(const char* name, size_t size)
	{
		int dim;

		checkError(nc_def_dim(identifier(), name, size, &dim));

		m_dimensions.push_back(Dimension(dim, name, size));

		return m_dimensions.back();
	}

	NetcdfEntity& createEntity(const char* name, const Type &type, size_t numDimensions, Dimension* dimensions)
	{
		m_entities.push_back(NetcdfEntity(name, type, m_ncDimSize, numDimensions, dimensions, offset(), *this));

		return m_entities.back();
	}

	/**
	 * @overload
	 */
	NetcdfEntity& createEntity(const char* name, const Type &type, std::vector<Dimension> &dimensions)
	{
		return createEntity(name, type, dimensions.size(), &dimensions[0]);
	}

	/**
	 * @overload
	 */
	NetcdfEntity& createEntity(const char* name, const Type &type)
	{
		return createEntity(name, type, 0, 0L);
	}

protected:
	bool setOffset(size_t partition)
	{

		if (partition >= numPartitions())
			// The last partition sets the offset for the first one
			partition = 0;

		unsigned long long o = offset()[partition];

		if (checkError(nc_put_var1_ulonglong(identifier(), m_ncVarOffset, &partition, &o)))
			return false;

		return true;
	}
};

}

#endif // PUML_NETCDF_GROUP_H
