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

#ifndef PUML_NETCDF_ENTITY_H
#define PUML_NETCDF_ENTITY_H

#include <vector>

#ifdef PARALLEL
#include <netcdf_par.h>
#endif // PARALLEL
#include <netcdf.h>

#include "PUML/Dimension.h"
#include "PUML/Entity.h"
#include "PUML/NetcdfElement.h"
#include "PUML/Type.h"

namespace PUML
{

class NetcdfEntity : public Entity, public NetcdfElement
{
	/** Partition dimension + user dimensions */
	std::vector<size_t> m_dimension;

public:
	NetcdfEntity()
	{
	}

	/**
	 * @param dimSize The netCDF dimension of group that contains the size
	 */
	NetcdfEntity(const char* name, const Type &type, int dimSize,
			size_t numUserDimensions, Dimension* userDimensions,
			const std::vector<size_t> &offset,
			NetcdfElement &group)
		: Entity(name, offset), NetcdfElement(&group),
		  m_dimension(numUserDimensions+1)
	{
		int ncVar;

		// Use std::vector to avoid memory leaks
		std::vector<int> dims(numUserDimensions+1);
		dims[0] = dimSize;
		for (size_t i = 0; i < numUserDimensions; i++) {
			dims[i+1] = userDimensions[i].identifier();

			// Set the size of the user dimension, we need them later
			m_dimension[i+1] = userDimensions[i].size();
		}

		if (checkError(nc_def_var(parentIdentifier(), name, type2nc(type), dims.size(), &dims[0], &ncVar)))
			return;

		setIdentifier(ncVar);
	}

#ifdef PARALLEL
	/**
	 * Overriding this function is only done in the parallel version
	 */
	bool setCollective(bool collective)
	{
		if (!Entity::setCollective(collective))
			return false;

		if (checkError(nc_var_par_access(parentIdentifier(), identifier(), (collective ? NC_COLLECTIVE : NC_INDEPENDENT))))
			return false;

		return true;
	}
#endif // PARALLEL

protected:
	bool _put(size_t partition, size_t size, const void* values)
	{
		// Use std::vector to avoid memory leaks
		std::vector<size_t> start(m_dimension.size(), 0);
		start[0] = offset(partition);

		// Set partition dimension for this call (not threadsafe)
		m_dimension[0] = size;

		if (checkError(nc_put_vara(parentIdentifier(), identifier(), &start[0], &m_dimension[0], values)))
			return false;

		return true;
	}

	bool _put_schar(size_t partition, size_t size, const signed char* values)
	{
		// Use std::vector to avoid memory leaks
		std::vector<size_t> start(m_dimension.size(), 0);
		start[0] = offset(partition);

		// Set partition dimension for this call (not threadsafe)
		m_dimension[0] = size;

		if (checkError(nc_put_vara_schar(parentIdentifier(), identifier(), &start[0], &m_dimension[0], values)))
			return false;

		return true;
	}

	bool _put_uchar(size_t partition, size_t size, const unsigned char* values)
	{
		// Use std::vector to avoid memory leaks
		std::vector<size_t> start(m_dimension.size(), 0);
		start[0] = offset(partition);

		// Set partition dimension for this call (not threadsafe)
		m_dimension[0] = size;

		if (checkError(nc_put_vara_uchar(parentIdentifier(), identifier(), &start[0], &m_dimension[0], values)))
			return false;

		return true;
	}

	bool _put_short(size_t partition, size_t size, const short* values)
	{
		// Use std::vector to avoid memory leaks
		std::vector<size_t> start(m_dimension.size(), 0);
		start[0] = offset(partition);

		// Set partition dimension for this call (not threadsafe)
		m_dimension[0] = size;

		if (checkError(nc_put_vara_short(parentIdentifier(), identifier(), &start[0], &m_dimension[0], values)))
			return false;

		return true;
	}

	bool _put_int(size_t partition, size_t size, const int* values)
	{
		// Use std::vector to avoid memory leaks
		std::vector<size_t> start(m_dimension.size(), 0);
		start[0] = offset(partition);

		// Set partition dimension for this call (not threadsafe)
		m_dimension[0] = size;

		if (checkError(nc_put_vara_int(parentIdentifier(), identifier(), &start[0], &m_dimension[0], values)))
			return false;

		return true;
	}

	bool _put_long(size_t partition, size_t size, const long* values)
	{
		// Use std::vector to avoid memory leaks
		std::vector<size_t> start(m_dimension.size(), 0);
		start[0] = offset(partition);

		// Set partition dimension for this call (not threadsafe)
		m_dimension[0] = size;

		if (checkError(nc_put_vara_long(parentIdentifier(), identifier(), &start[0], &m_dimension[0], values)))
			return false;

		return true;
	}

	bool _put_float(size_t partition, size_t size, const float* values)
	{
		// Use std::vector to avoid memory leaks
		std::vector<size_t> start(m_dimension.size(), 0);
		start[0] = offset(partition);

		// Set partition dimension for this call (not threadsafe)
		m_dimension[0] = size;

		if (checkError(nc_put_vara_float(parentIdentifier(), identifier(), &start[0], &m_dimension[0], values)))
			return false;

		return true;
	}

	bool _put_double(size_t partition, size_t size, const double* values)
	{
		// Use std::vector to avoid memory leaks
		std::vector<size_t> start(m_dimension.size(), 0);
		start[0] = offset(partition);

		// Set partition dimension for this call (not threadsafe)
		m_dimension[0] = size;

		if (checkError(nc_put_vara_double(parentIdentifier(), identifier(), &start[0], &m_dimension[0], values)))
			return false;

		return true;
	}

	bool _put_ushort(size_t partition, size_t size, const unsigned short* values)
	{
		// Use std::vector to avoid memory leaks
		std::vector<size_t> start(m_dimension.size(), 0);
		start[0] = offset(partition);

		// Set partition dimension for this call (not threadsafe)
		m_dimension[0] = size;

		if (checkError(nc_put_vara_ushort(parentIdentifier(), identifier(), &start[0], &m_dimension[0], values)))
			return false;

		return true;
	}

	bool _put_uint(size_t partition, size_t size, const unsigned int* values)
	{
		// Use std::vector to avoid memory leaks
		std::vector<size_t> start(m_dimension.size(), 0);
		start[0] = offset(partition);

		// Set partition dimension for this call (not threadsafe)
		m_dimension[0] = size;

		if (checkError(nc_put_vara_uint(parentIdentifier(), identifier(), &start[0], &m_dimension[0], values)))
			return false;

		return true;
	}

	bool _put_longlong(size_t partition, size_t size, const long long* values)
	{
		// Use std::vector to avoid memory leaks
		std::vector<size_t> start(m_dimension.size(), 0);
		start[0] = offset(partition);

		// Set partition dimension for this call (not threadsafe)
		m_dimension[0] = size;

		if (checkError(nc_put_vara_longlong(parentIdentifier(), identifier(), &start[0], &m_dimension[0], values)))
			return false;

		return true;
	}

	bool _put_ulonglong(size_t partition, size_t size, const unsigned long long* values)
	{
		// Use std::vector to avoid memory leaks
		std::vector<size_t> start(m_dimension.size(), 0);
		start[0] = offset(partition);

		// Set partition dimension for this call (not threadsafe)
		m_dimension[0] = size;

		if (checkError(nc_put_vara_ulonglong(parentIdentifier(), identifier(), &start[0], &m_dimension[0], values)))
			return false;

		return true;
	}

private:
	/**
	 * @return The nc identifier for this type
	 */
	static nc_type type2nc(const Type &type)
	{
		switch (type.baseType()) {
		case Type::CHAR:
			return NC_CHAR;
		case Type::BYTE:
			return NC_BYTE;
		case Type::SHORT:
			return NC_SHORT;
		case Type::INT:
			return NC_INT;
		case Type::INT64:
			return NC_INT64;
		case Type::FLOAT:
			return NC_FLOAT;
		case Type::DOUBLE:
			return NC_DOUBLE;
		case Type::UBYTE:
			return NC_UBYTE;
		case Type::UINT:
			return NC_UINT;
		case Type::UINT64:
			return NC_UINT64;
		default:
			return type.identifier();
		}
	}
};

}

#endif // PUML_ENTITY_H
