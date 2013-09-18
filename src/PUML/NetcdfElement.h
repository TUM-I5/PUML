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

#ifndef PUML_NETCDF_ELEMENT_H
#define PUML_NETCDF_ELEMENT_H

#include <string>

#include <netcdf.h>

namespace PUML
{

/**
 * Contains basic functionality to handle a nc element (file, group, variable, ...)
 */
class NetcdfElement
{
private:
	/** nc identifier */
	int m_ncIdentifier;

	/** Parent element */
	NetcdfElement* m_parent;

	/** nc error (or NC_NOERR if no error occurred) */
	int m_ncError;

public:
	NetcdfElement(NetcdfElement* parent = 0L)
		: m_ncIdentifier(-1), m_parent(parent), m_ncError(NC_NOERR)
	{
	}

	NetcdfElement(int identifier, NetcdfElement* parent = 0L)
		: m_ncIdentifier(identifier), m_parent(parent), m_ncError(NC_NOERR)
	{
	}

	virtual ~NetcdfElement()
	{
	}

	/**
	 * @return The netCDF identifier
	 *
	 * @ingroup LowLevelApi
	 */
	int identifier() const
	{
		return m_ncIdentifier;
	}

	/**
	 * @return True if no error occurred, false otherwise
	 */
	bool isValid() const
	{
		if (m_parent)
			return m_parent->isValid();

		return m_ncError == NC_NOERR;
	}

	/**
	 * @return The message for the error
	 */
	std::string errorMsg() const
	{
		if (m_parent)
			return m_parent->errorMsg();

		return nc_strerror(m_ncError);
	}

protected:
	void setIdentifier(int identifier)
	{
		m_ncIdentifier = identifier;
	}

	/**
	 * Returns the identifier of the parent if a parent exists
	 *
	 * @return The netCDF identifier of the parent element
	 */
	int parentIdentifier() const
	{
		if (m_parent)
			return m_parent->m_ncIdentifier;

		return m_ncIdentifier;
	}

	/**
	 * @return the nc identifier of the file (top most identifier)
	 */
	int fileIdentifier() const
	{
		if (m_parent)
			return m_parent->fileIdentifier();

		return m_ncIdentifier;
	}

	/**
	 * Checks if result contains an error and saves the error state
	 *
	 * @return True result is an error false otherwise
	 */
	bool checkError(int result)
	{
		if (m_parent)
			// Propagate error to the parent element
			return m_parent->checkError(result);

		if (m_ncError == NC_NOERR)
			// An error occurred early -> ignore this error
			m_ncError = result;

		return result != NC_NOERR;
	}
};

}

#endif // PUML_NETCDF_ELEMENT_H
