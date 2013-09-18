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

#ifndef PUML_DIMENSION_H
#define PUML_DIMENSION_H

#include <string>

namespace PUML
{

class Dimension
{
private:
	/** Identifier (for internal use) */
	long m_identifier;

	/** Name of the dimension */
	std::string m_name;

	/** Size of the dimension */
	size_t m_size;

public:
	Dimension(long identifier, const char* name, size_t size)
		: m_identifier(identifier), m_name(name), m_size(size)
	{
	}

	long identifier() const
	{
		return m_identifier;
	}

	size_t size() const
	{
		return m_size;
	}
};

}

#endif // PUML_DIMENSION_H
