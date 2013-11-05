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
public:
	NetcdfEntity()
	{
	}

	/**
	 * @param dimSize The netCDF dimension of the group that contains the size
	 */
	NetcdfEntity(const char* name, const Type &type, int dimSize,
			size_t numUserDimensions, const Dimension* userDimensions,
			const std::vector<size_t> &offset, NetcdfEntity* index,
			NetcdfElement &group, MPIElement &comm)
		: Entity(name, numUserDimensions, userDimensions, offset, index, comm), NetcdfElement(&group)
	{
		int ncVar;

		// Use std::vector to avoid memory leaks
		std::vector<int> dims(numUserDimensions+1);
		dims[0] = dimSize;
		for (size_t i = 0; i < numUserDimensions; i++)
			dims[i+1] = userDimensions[i].identifier();

		if (checkError(nc_def_var(parentIdentifier(), name, type2nc(type), dims.size(), &dims[0], &ncVar)))
			return;

		setIdentifier(ncVar);
	}

	/**
	 * Constructor to load an entity from a nc file
	 */
	NetcdfEntity(int ncId, const std::vector<size_t> &offset, NetcdfEntity* index, NetcdfElement &group, MPIElement &comm)
		: Entity(offset, index, comm), NetcdfElement(ncId, &group)
	{
		char name[NC_MAX_NAME+1];
		if (checkError(nc_inq_varname(parentIdentifier(), identifier(), name)))
			return;
		setName(name);

		// Learn about the dimension size
		int numDims;
		if (checkError(nc_inq_varndims(parentIdentifier(), identifier(), &numDims)))
			return;

		std::vector<int> dims(numDims);
		if (checkError(nc_inq_vardimid(parentIdentifier(), identifier(), &dims[0])))
			return;

		dimSize().resize(numDims);
		for (int i = 1; i < numDims; i++) { // Skip pum dimension
			if (checkError(nc_inq_dimlen(parentIdentifier(), dims[i], &dimSize()[i])))
				return;
		}
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
	bool _puta(const size_t* start, const size_t* size, const void* values)
	{
		return !checkError(nc_put_vara(parentIdentifier(), identifier(), start, size, values));
	}

	bool _puta_schar(const size_t* start, const size_t* size, const signed char* values)
	{
		return !checkError(nc_put_vara_schar(parentIdentifier(), identifier(), start, size, values));
	}

	bool _puta_uchar(const size_t* start, const size_t* size, const unsigned char* values)
	{
		return !checkError(nc_put_vara_uchar(parentIdentifier(), identifier(), start, size, values));
	}

	bool _puta_short(const size_t* start, const size_t* size, const short* values)
	{
		return !checkError(nc_put_vara_short(parentIdentifier(), identifier(), start, size, values));
	}

	bool _puta_int(const size_t* start, const size_t* size, const int* values)
	{
		return !checkError(nc_put_vara_int(parentIdentifier(), identifier(), start, size, values));
	}

	bool _puta_long(const size_t* start, const size_t* size, const long* values)
	{
		return !checkError(nc_put_vara_long(parentIdentifier(), identifier(), start, size, values));
	}

	bool _puta_float(const size_t* start, const size_t* size, const float* values)
	{
		return !checkError(nc_put_vara_float(parentIdentifier(), identifier(), start, size, values));
	}

	bool _puta_double(const size_t* start, const size_t* size, const double* values)
	{
		return !checkError(nc_put_vara_double(parentIdentifier(), identifier(), start, size, values));
	}

	bool _puta_ushort(const size_t* start, const size_t* size, const unsigned short* values)
	{
		return !checkError(nc_put_vara_ushort(parentIdentifier(), identifier(), start, size, values));
	}

	bool _puta_uint(const size_t* start, const size_t* size, const unsigned int* values)
	{
		return !checkError(nc_put_vara_uint(parentIdentifier(), identifier(), start, size, values));
	}

	bool _puta_longlong(const size_t* start, const size_t* size, const long long* values)
	{
		return !checkError(nc_put_vara_longlong(parentIdentifier(), identifier(), start, size, values));
	}

	bool _puta_ulonglong(const size_t* start, const size_t* size, const unsigned long long* values)
	{
		return !checkError(nc_put_vara_ulonglong(parentIdentifier(), identifier(), start, size, values));
	}

	bool _geta(const size_t* start, const size_t* size, void* values)
	{
		return !checkError(nc_get_vara(parentIdentifier(), identifier(), start, size, values));
	}

	bool _geta_schar(const size_t* start, const size_t* size, signed char* values)
	{
		return !checkError(nc_get_vara_schar(parentIdentifier(), identifier(), start, size, values));
	}

	bool _geta_uchar(const size_t* start, const size_t* size, unsigned char* values)
	{
		return !checkError(nc_get_vara_uchar(parentIdentifier(), identifier(), start, size, values));
	}

	bool _geta_short(const size_t* start, const size_t* size, short* values)
	{
		return !checkError(nc_get_vara_short(parentIdentifier(), identifier(), start, size, values));
	}

	bool _geta_int(const size_t* start, const size_t* size, int* values)
	{
		return !checkError(nc_get_vara_int(parentIdentifier(), identifier(), start, size, values));
	}

	bool _geta_long(const size_t* start, const size_t* size, long* values)
	{
		return !checkError(nc_get_vara_long(parentIdentifier(), identifier(), start, size, values));
	}

	bool _geta_float(const size_t* start, const size_t* size, float* values)
	{
		return !checkError(nc_get_vara_float(parentIdentifier(), identifier(), start, size, values));
	}

	bool _geta_double(const size_t* start, const size_t* size, double* values)
	{
		return !checkError(nc_get_vara_double(parentIdentifier(), identifier(), start, size, values));
	}

	bool _geta_ushort(const size_t* start, const size_t* size, unsigned short* values)
	{
		return !checkError(nc_get_vara_ushort(parentIdentifier(), identifier(), start, size, values));
	}

	bool _geta_uint(const size_t* start, const size_t* size, unsigned int* values)
	{
		return !checkError(nc_get_vara_uint(parentIdentifier(), identifier(), start, size, values));
	}

	bool _geta_longlong(const size_t* start, const size_t* size, long long* values)
	{
		return !checkError(nc_get_vara_longlong(parentIdentifier(), identifier(), start, size, values));
	}

	bool _geta_ulonglong(const size_t* start, const size_t* size, unsigned long long* values)
	{
		return !checkError(nc_get_vara_ulonglong(parentIdentifier(), identifier(), start, size, values));
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
