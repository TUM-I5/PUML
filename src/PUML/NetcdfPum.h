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

	NetcdfGroup& createGroup(const char* name, size_t size = Group::UNLIMITED)
	{
		if (size == Group::UNLIMITED)
			size = NC_UNLIMITED;

		m_groups.push_back(NetcdfGroup(name, numPartitions(), *this, *this, size));

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
};

}

#endif // PUML_NETCDF_PUM_H
