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

#ifndef PUML_NETCDF_PUM_H
#define PUML_NETCDF_PUM_H

#ifdef PARALLEL
#include <mpi.h>
#endif // PARALLEL

#include <vector>

#ifdef PARALLEL
#include <netcdf_par.h>
//#endif // PARALLEL
#else
#include <netcdf.h>
#endif

#include "PUML/NetcdfGroup.h"
#include "PUML/NetcdfElement.h"
#include "PUML/Pum.h"

namespace PUML
{

class NetcdfPum : public Pum, public NetcdfElement
{
private:
	/** Groups in this file */
	std::vector<NetcdfGroup> m_groups;

public:
	NetcdfPum()
	{
	}

	virtual ~NetcdfPum()
	{
		nc_close(identifier());
	}

	bool open(const char* path)
	{
		int ncFile;

		if (checkError(nc_open(path, NC_NETCDF4, &ncFile)))
				return false;
		setIdentifier(ncFile);

		return loadFile();
	}

#ifdef PARALLEL
	/**
	 * Overridden to work around overload/subclass issues
	 */
	bool open(const char* path, MPI_Comm comm, MPI_Info info = MPI_INFO_NULL)
	{
		return Pum::open(path, comm, info);
	}
#endif // PARALLEL

	NetcdfGroup& createGroup(const char* name, size_t size = Group::UNLIMITED)
	{
		if (size == Group::UNLIMITED)
			size = NC_UNLIMITED;

		m_groups.push_back(NetcdfGroup(name, numPartitions(), *this, *this, size));

		return m_groups.back();
	}

	NetcdfGroup& createGroupIndexed(const char* name, size_t size = Group::UNLIMITED, size_t indexSize = Group::UNLIMITED)
	{
		if (size == Group::UNLIMITED)
			size = NC_UNLIMITED;
		if (indexSize == Group::UNLIMITED)
			indexSize = NC_UNLIMITED;

		m_groups.push_back(NetcdfGroup(name, numPartitions(), *this, *this, size, indexSize));

		return m_groups.back();
	}

	bool endDefinition()
	{
		if (!Pum::endDefinition())
			return false;

		return !checkError(nc_enddef(identifier()));
	}

	bool close()
	{
		if (checkError(nc_close(identifier())))
			return false;

		return true;
	}

protected:
	bool _create(const char* path)
	{
		int ncFile;

		if (checkError(nc_create(path, NC_NETCDF4, &ncFile)))
			return false;
		setIdentifier(ncFile);

		return initFile();
	}

#ifdef PARALLEL
	bool _create(const char* path, MPI_Comm comm, MPI_Info info = MPI_INFO_NULL)
	{
		int ncFile;

		if (checkError(nc_create_par(path, NC_NETCDF4 | NC_MPIIO, comm, info, &ncFile)))
				return false;
		setIdentifier(ncFile);

		return initFile();
	}

	bool _open(const char* path, MPI_Comm comm, MPI_Info info = MPI_INFO_NULL)
	{
		int ncFile;

		if (checkError(nc_open_par(path, NC_NETCDF4 | NC_MPIIO, comm, info, &ncFile)))
				return false;
		setIdentifier(ncFile);

		return loadFile();
	}
#endif // PARALLEL

private:
	/**
	 * Initialize a new nc pum file
	 */
	bool initFile()
	{
		if (checkError(nc_put_att_text(identifier(), NC_GLOBAL, ATT_CONVENTIONS, CONVENTIONS.size(), CONVENTIONS.c_str())))
			return false;
		if (checkError(nc_put_att_int(identifier(), NC_GLOBAL, ATT_FILE_VERSION, NC_INT, 1, &FILE_VERSION)))
			return false;

		return true;
	}

	/**
	 * Check nc pum file and load groups, etc
	 */
	bool loadFile()
	{
		size_t len;
		if (checkError(nc_inq_attlen(identifier(), NC_GLOBAL, ATT_CONVENTIONS, &len)))
			return false;

		std::vector<char> conventions(len); // Use std::vector to avoid memory leaks
		if (checkError(nc_get_att_text(identifier(), NC_GLOBAL, ATT_CONVENTIONS, &conventions[0])))
			return false;
		if (CONVENTIONS.compare(&conventions[0]) != 0)
			return false;

		int fileVersion;
		if (checkError(nc_get_att_int(identifier(), NC_GLOBAL, ATT_FILE_VERSION, &fileVersion)))
			return false;
		if (fileVersion != FILE_VERSION)
			// Currently only one version is supported
			return false;

		// Get group ids
		int numGrps;
		if (checkError(nc_inq_grps(identifier(), &numGrps, 0L)))
			return false;
		std::vector<int> grpIds(numGrps);
		if (checkError(nc_inq_grps(identifier(), 0L, &grpIds[0])))
			return false;

		// Create the groups
		for (std::vector<int>::const_iterator i = grpIds.begin(); i != grpIds.end(); i++)
			m_groups.push_back(NetcdfGroup(*i, *this, *this));

		return true;

	}
};

}

#endif // PUML_NETCDF_PUM_H
