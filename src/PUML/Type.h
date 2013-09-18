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

#ifndef PUML_TYPE_H
#define PUML_TYPE_H

namespace PUML
{

class Type
{
public:
	/**
	 * Base types for entities (a primitive or custom type)
	 */
	enum BaseType
	{
		CHAR,
		BYTE,
		SHORT,
		INT,
		INT64,
		FLOAT,
		DOUBLE,
		UBYTE,
		USHORT,
		UINT,
		UINT64,
		CUSTOM
	};

private:
	/** Primitive or CUSTOM */
	BaseType m_type;

	/** Identifier for a custom type */
	long m_customIdentifier;

public:
	/**
	 * Create a primitive type
	 */
	Type(BaseType type)
		: m_type(type), m_customIdentifier(-1)
	{
	}

	/**
	 * Create a custom type
	 *
	 * @ingroup LowLevelApi
	 */
	Type(long identifier)
		: m_type(CUSTOM), m_customIdentifier(identifier)
	{
	}

	BaseType baseType() const
	{
		return m_type;
	}

	long identifier() const
	{
		return m_customIdentifier;
	}

public:
	static const Type Char;
	static const Type Byte;
	static const Type Short;
	static const Type Int;
	static const Type Int64;
	static const Type Float;
	static const Type Double;
	static const Type uByte;
	static const Type uShort;
	static const Type uInt;
	static const Type uInt64;
};

}

#endif // PUML_TYPE_H
