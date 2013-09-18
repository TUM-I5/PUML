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

#ifndef PUML_ENTITY_H
#define PUML_ENTITY_H

#include <limits>
#include <string>

namespace PUML
{

class Entity
{
private:
	/** Name of this entity */
	std::string m_name;

	/** This entity currently in collective mode? */
	bool m_collective;

	/** Reference to the offset of the group */
	const std::vector<size_t>* m_offset;

public:
	Entity()
		: m_collective(false), m_offset(0L)
	{
	}

	Entity(const char* name, const std::vector<size_t> &offset)
		: m_name(name), m_collective(false), m_offset(&offset)
	{
	}

	virtual ~Entity()
	{
	}

	/**
	 * Activate/Deactivate collective mode
	 */
	virtual bool setCollective(bool collective)
	{
		m_collective = collective;
		return true;
	}

	/**
	 * Writes the values for one partition the file. Make sure that the size of the partition is already determined.
	 *
	 * @param size Number of elements that should be stored
	 *
	 * @see Group::setSize
	 */
	template<typename T>
	bool put(size_t partition, size_t size, const T* values)
	{
		if (!isPartitionSizeSet(partition))
			return false;

		return __put(partition, size, values);
	}

protected:
	size_t offset(size_t partition)
	{
		return (*m_offset)[partition];
	}

	virtual bool _put(size_t partition, size_t size, const void* values) = 0;
	virtual bool _put_schar(size_t partition, size_t size, const signed char* values) = 0;
	virtual bool _put_uchar(size_t partition, size_t size, const unsigned char* values) = 0;
	virtual bool _put_short(size_t partition, size_t size, const short* values) = 0;
	virtual bool _put_int(size_t partition, size_t size, const int* values) = 0;
	virtual bool _put_long(size_t partition, size_t size, const long* values) = 0;
	virtual bool _put_float(size_t partition, size_t size, const float* values) = 0;
	virtual bool _put_double(size_t partition, size_t size, const double* values) = 0;
	virtual bool _put_ushort(size_t partition, size_t size, const unsigned short* values) = 0;
	virtual bool _put_uint(size_t partition, size_t size, const unsigned int* values) = 0;
	virtual bool _put_longlong(size_t partition, size_t size, const long long* values) = 0;
	virtual bool _put_ulonglong(size_t partition, size_t size, const unsigned long long* values) = 0;

private:
	bool isPartitionSizeSet(size_t partition)
	{
		return (*m_offset)[partition] != std::numeric_limits<size_t>::max();
	}

	template<typename T>
	bool __put(size_t partition, size_t size, const T* values)
	{
		return _put(partition, size, values);
	}
};

template<> inline
bool Entity::__put(size_t partition, size_t size, const signed char* values)
{ return _put_schar(partition, size, values); }

template<> inline
bool Entity::__put(size_t partition, size_t size, const unsigned char* values)
{ return _put_uchar(partition, size, values); }

template<> inline
bool Entity::__put(size_t partition, size_t size, const short* values)
{ return _put_short(partition, size, values); }

template<> inline
bool Entity::__put(size_t partition, size_t size, const int* values)
{ return _put_int(partition, size, values); }

template<> inline
bool Entity::__put(size_t partition, size_t size, const long* values)
{ return _put_long(partition, size, values); }

template<> inline
bool Entity::__put(size_t partition, size_t size, const float* values)
{ return _put_float(partition, size, values); }

template<> inline
bool Entity::__put(size_t partition, size_t size, const double* values)
{ return _put_double(partition, size, values); }

template<> inline
bool Entity::__put(size_t partition, size_t size, const unsigned short* values)
{ return _put_ushort(partition, size, values); }

template<> inline
bool Entity::__put(size_t partition, size_t size, const unsigned int* values)
{ return _put_uint(partition, size, values); }

template<> inline
bool Entity::__put(size_t partition, size_t size, const long long* values)
{ return _put_longlong(partition, size, values); }

template<> inline
bool Entity::__put(size_t partition, size_t size, const unsigned long long* values)
{ return _put_ulonglong(partition, size, values); }

}

#endif // PUML_ENTITY_H
