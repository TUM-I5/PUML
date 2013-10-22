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

#include <map>

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

	/** nc identifier for the index size dimension */
	int m_ncDimIndexSize;

	/** User defined dimensions */
	std::vector<Dimension> m_dimensions;

	/** nc identifier for the offset variable */
	int m_ncVarOffset;

	/** index variable */
	NetcdfEntity m_entityIndex;

	/** Entities in this group */
	std::map<std::string, NetcdfEntity> m_entities;

public:
	NetcdfGroup()
		: m_ncDimPartition(-1), m_ncDimSize(-1), m_ncDimIndexSize(-1), m_ncVarOffset(-1)
	{
	}

	/**
	 * @param size The total size of this group. Use NC_UNLIMITED if unknown
	 */
	NetcdfGroup(const char* name, size_t numPartitions, NetcdfElement &ncPum, MPIElement &comm, size_t size)
		: Group(name, numPartitions, comm), NetcdfElement(&ncPum), m_ncDimIndexSize(-1)
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

	/**
	 * Constructor to load a group from the nc file
	 *
	 * @param ncId The netcdf identifier of the group
	 */
	NetcdfGroup(int ncId, NetcdfElement &ncPum, MPIElement &comm)
		: Group(comm), NetcdfElement(ncId, &ncPum)
	{
		char name[NC_MAX_NAME+1];
		if (checkError(nc_inq_grpname(identifier(), name)))
			return;
		setName(name);

		if (checkError(nc_inq_dimid(identifier(), DIM_PARTITION, &m_ncDimPartition)))
			return;

		if (checkError(nc_inq_dimid(identifier(), DIM_SIZE, &m_ncDimSize)))
			return;

		int ncError = nc_inq_dimid(identifier(), DIM_INDEXSIZE, &m_ncDimIndexSize);
		if (ncError == NC_EBADDIM) {
			// Not an indexed variable
			m_ncDimIndexSize = -1;
		} else {
			if (checkError(ncError))
				return;

			// TODO Indexed group
		}

		if (checkError(nc_inq_varid(identifier(), VAR_OFFSET, &m_ncVarOffset)))
			return;
#ifdef PARALLEL
		if (checkError(nc_var_par_access(identifier(), m_ncVarOffset, NC_COLLECTIVE)))
			return;
#endif // PARALLEL

		// Read offsets
		size_t numPartitions;
		if (checkError(nc_inq_dimlen(identifier(), m_ncDimPartition, &numPartitions)))
			return;

		std::vector<unsigned long long> o(numPartitions);
		if (checkError(nc_get_var_ulonglong(identifier(), m_ncVarOffset, &o[0])))
			return;

		offset().resize(numPartitions+1);
		for (size_t i = 0; i < numPartitions; i++)
			offset()[i] = o[i];

		if (checkError(nc_inq_dimlen(identifier(), m_ncDimSize, &offset()[numPartitions])))
			return;
	}

	Dimension& createDimension(const char* name, size_t size)
	{
		int dim;

		checkError(nc_def_dim(identifier(), name, size, &dim));

		m_dimensions.push_back(Dimension(dim, name, size));

		return m_dimensions.back();
	}

	NetcdfEntity* createEntity(const char* name, const Type &type, size_t numDimensions, Dimension* dimensions)
	{
		NetcdfEntity entity = NetcdfEntity(name, type, m_ncDimSize, numDimensions, dimensions, offset(), *this);
		if (!entity.isValid())
			return 0L;

		m_entities[name] = entity;

		return &m_entities[name];
	}

	/**
	 * @overload
	 */
	NetcdfEntity* createEntity(const char* name, const Type &type, std::vector<Dimension> &dimensions)
	{
		return createEntity(name, type, dimensions.size(), &dimensions[0]);
	}

	/**
	 * @overload
	 */
	NetcdfEntity* createEntity(const char* name, const Type &type)
	{
		return createEntity(name, type, 0, 0L);
	}

	NetcdfEntity* getEntity(const char* name)
	{
		if (m_entities.find(name) == m_entities.end())
			return 0L;

		return &m_entities.at(name);
	}

	/**
	 * Loads the entities from the netcdf file
	 * We can't do this in the constructor because this results in wrong values for m_parent
	 *
	 * @internal
	 */
	bool loadEntities()
	{
		int numVars;
		if (checkError(nc_inq_varids(identifier(), &numVars, 0L)))
			return false;

		std::vector<int> varIds(numVars);
		if (checkError(nc_inq_varids(identifier(), 0L, &varIds[0])))
			return false;

		for (std::vector<int>::const_iterator i = varIds.begin(); i != varIds.end(); i++) {
			if (*i == m_ncVarOffset)
				continue;

			NetcdfEntity entity = NetcdfEntity(*i, offset(), *this);
			m_entities[entity.name()] = entity;
		}

		// TODO index entity

		return true;
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

	NetcdfEntity* _addIndex(size_t indexSize)
	{
		if (checkError(nc_def_dim(identifier(), DIM_INDEXSIZE, indexSize, &m_ncDimIndexSize)))
			return 0L;

		m_entityIndex = NetcdfEntity(VAR_INDEX, Type::UINT64, m_ncDimIndexSize, 0, 0L, offset(), *this);

		return &m_entityIndex;
	}
};

}

#endif // PUML_NETCDF_GROUP_H
