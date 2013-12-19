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

#include "PUML/MPIElement.h"

namespace PUML
{

class Entity : protected MPIElement
{
private:
	/** Name of this entity */
	std::string m_name;

	/** This entity currently in collective mode? */
	bool m_collective;

	/** Size of partition dimension (temporary) + user dimensions */
	std::vector<size_t> m_dimSize;

	/** Pointer to the offset of the group */
	const std::vector<size_t>* m_offset;

	/** Pointer to the index entity */
	Entity* m_index;

	/**
	 * Helper structure to sum up contiguous indices
	 */
	struct IndexedRange
	{
		/** Position in the netCDF file */
		size_t pos;
		/** Number of continues values */
		size_t count;
		/** Position in the local copy of the values */
		size_t localPos;
	};

public:
	Entity()
		: m_collective(false), m_offset(0L), m_index(0L)
	{
	}

	Entity(const char* name, size_t numUserDimensions, const Dimension* userDimensions,
			const std::vector<size_t> &offset, Entity* index, MPIElement &comm)
		: MPIElement(comm),
		  m_name(name), m_collective(false),
		  m_dimSize(numUserDimensions+1), m_offset(&offset), m_index(index)
	{
		for (size_t i = 0; i < numUserDimensions; i++) {
			// Set the size of the user dimension, we need them later
			m_dimSize[i+1] = userDimensions[i].size();
		}
	}

	/**
	 * Constructor for loading an entity from file
	 */
	Entity(const std::vector<size_t> &offset, Entity* index, MPIElement &comm)
		: MPIElement(comm),
		  m_collective(false), m_offset(&offset), m_index(index)
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
		if (!isPartitionOffsetSet(partition))
			return false;

		if (!indexed())
			return puta((*m_offset)[partition], size, values);

		// compute position and count of values
		std::vector<IndexedRange> valuePos;
		size_t accesses;
		if (!getValuePos(partition, size, valuePos, accesses))
			return false;


		for (size_t i = 0; i < accesses; i++) {
			// Due to collective I/O accesses might be larger than valuePos.size()
			IndexedRange& v = valuePos[i % valuePos.size()];
			if (!puta(v.pos, v.count, &values[v.localPos]))
				return false;
		}

		return true;
	}

	template<typename T>
	bool get(size_t partition, size_t size, T* values)
	{
		if (!isPartitionOffsetSet(partition))
			return false;

		if (!indexed())
			return geta((*m_offset)[partition], size, values);

		// compute position and count of values
		std::vector<IndexedRange> valuePos;
		size_t accesses;
		if (!getValuePos(partition, size, valuePos, accesses))
			return false;

		for (size_t i = 0; i < accesses; i++) {
			// Due to collective I/O accesses might be larger than valuePos.size()
			IndexedRange& v = valuePos[i % valuePos.size()];
			if (!geta(v.pos, v.count, &values[v.localPos]))
				return false;
		}

		return true;
	}

	/**
	 * @return All values of a partition
	 */
	template<typename T>
	bool get(size_t partition, T* values)
	{
		if (!isPartitionSizeSet(partition))
			return false;

		return get(partition, partitionSize(partition), values);
	}

	/**
	 * Put values at absolute position
	 */
	template<typename T>
	bool puta(size_t start, size_t size, const T* values)
	{
		// Use std::vector to avoid memory leaks
		std::vector<size_t> s(m_dimSize.size(), 0);
		s[0] = start;

		// Set partition dimension for this call (not threadsafe)
		m_dimSize[0] = size;

		return __puta(&s[0], &m_dimSize[0], values);
	}

	/**
	 * Get values at absolute position
	 */
	template<typename T>
	bool geta(size_t start, size_t size, T* values)
	{
		// Use std::vector to avoid memory leaks
		std::vector<size_t> s(m_dimSize.size(), 0);
		s[0] = start;

		// Set partition dimension for this call (not threadsafe)
		m_dimSize[0] = size;

		return _geta(&s[0], &m_dimSize[0], values);
	}

	const char* name() const
	{
		return m_name.c_str();
	}

protected:
	/**
	 * @return The number of dimension of this entity
	 */
	int numDims() const
	{
		return m_dimSize.size();
	}

	const std::vector<size_t>& dimSize() const
	{
		return m_dimSize;
	}

	std::vector<size_t>& dimSize()
	{
		return m_dimSize;
	}

	void setName(const char* name)
	{
		m_name = name;
	}

	virtual bool _puta(const size_t* start, const size_t* size, const void* values) = 0;
	virtual bool _puta_schar(const size_t* start, const size_t* size, const signed char* values) = 0;
	virtual bool _puta_uchar(const size_t* start, const size_t* size, const unsigned char* values) = 0;
	virtual bool _puta_short(const size_t* start, const size_t* size, const short* values) = 0;
	virtual bool _puta_int(const size_t* start, const size_t* size, const int* values) = 0;
	virtual bool _puta_long(const size_t* start, const size_t* size, const long* values) = 0;
	virtual bool _puta_float(const size_t* start, const size_t* size, const float* values) = 0;
	virtual bool _puta_double(const size_t* start, const size_t* size, const double* values) = 0;
	virtual bool _puta_ushort(const size_t* start, const size_t* size, const unsigned short* values) = 0;
	virtual bool _puta_uint(const size_t* start, const size_t* size, const unsigned int* values) = 0;
	virtual bool _puta_longlong(const size_t* start, const size_t* size, const long long* values) = 0;
	virtual bool _puta_ulonglong(const size_t* start, const size_t* size, const unsigned long long* values) = 0;

	virtual bool _geta(const size_t* start, const size_t* size, void* values) = 0;
	virtual bool _geta_schar(const size_t* start, const size_t* size, signed char* values) = 0;
	virtual bool _geta_uchar(const size_t* start, const size_t* size, unsigned char* values) = 0;
	virtual bool _geta_short(const size_t* start, const size_t* size, short* values) = 0;
	virtual bool _geta_int(const size_t* start, const size_t* size, int* values) = 0;
	virtual bool _geta_long(const size_t* start, const size_t* size, long* values) = 0;
	virtual bool _geta_float(const size_t* start, const size_t* size, float* values) = 0;
	virtual bool _geta_double(const size_t* start, const size_t* size, double* values) = 0;
	virtual bool _geta_ushort(const size_t* start, const size_t* size, unsigned short* values) = 0;
	virtual bool _geta_uint(const size_t* start, const size_t* size, unsigned int* values) = 0;
	virtual bool _geta_longlong(const size_t* start, const size_t* size, long long* values) = 0;
	virtual bool _geta_ulonglong(const size_t* start, const size_t* size, unsigned long long* values) = 0;

private:
	bool isPartitionOffsetSet(size_t partition)
	{
		return (*m_offset)[partition] != std::numeric_limits<size_t>::max();
	}

	bool isPartitionSizeSet(size_t partition)
	{
		return (*m_offset)[partition+1] != std::numeric_limits<size_t>::max();
	}

	/**
	 * @return The size of a partition (if set)
	 *
	 * @see Group::size
	 */
	size_t partitionSize(size_t partition)
	{
		return (*m_offset)[partition+1] - (*m_offset)[partition];
	}

	bool indexed() const
	{
		return m_index != 0L;
	}

	template<typename T>
	bool __puta(const size_t* start, const size_t* size, const T* values)
	{
		return _puta(start, size, values);
	}

	template<typename T>
	bool __geta(const size_t* start, const size_t* size, T* values)
	{
		return _geta(start, size, values);
	}

	/**
	 * @param accesses The number of accesses we need (may differ from valuePos.size())
	 */
	bool getValuePos(size_t partition, size_t size, std::vector<IndexedRange> &valuePos, size_t &accesses)
	{
		std::vector<unsigned long> index(size);
		if (!m_index->get(partition, size, &index[0]))
			return false;

		IndexedRange range = {index[0], 1, 0};
		valuePos.push_back(range);
		for (size_t i = 1; i < index.size(); i++) {
			if (index[i] == valuePos.back().pos + valuePos.back().count)
				valuePos.back().count++;
			else {
				size_t localPos = valuePos.back().localPos + valuePos.back().count;
				IndexedRange range2 = {index[i], 1, localPos};
				valuePos.push_back(range2);
			}
		}

		// Get maximum number of access we need to get all data
		unsigned long maxAccesses = valuePos.size();
#ifdef PARALLEL
		if (m_collective)
			MPI_Allreduce(MPI_IN_PLACE, &maxAccesses, 1, MPI_UNSIGNED_LONG, MPI_MAX, mpiComm());
#endif // PARALLEL

		accesses = maxAccesses;

		return true;
	}
};

template<> inline
bool Entity::__puta(const size_t* start, const size_t* size, const signed char* values)
{ return _puta_schar(start, size, values); }

template<> inline
bool Entity::__puta(const size_t* start, const size_t* size, const unsigned char* values)
{ return _puta_uchar(start, size, values); }

template<> inline
bool Entity::__puta(const size_t* start, const size_t* size, const short* values)
{ return _puta_short(start, size, values); }

template<> inline
bool Entity::__puta(const size_t* start, const size_t* size, const int* values)
{ return _puta_int(start, size, values); }

template<> inline
bool Entity::__puta(const size_t* start, const size_t* size, const long* values)
{ return _puta_long(start, size, values); }

template<> inline
bool Entity::__puta(const size_t* start, const size_t* size, const float* values)
{ return _puta_float(start, size, values); }

template<> inline
bool Entity::__puta(const size_t* start, const size_t* size, const double* values)
{ return _puta_double(start, size, values); }

template<> inline
bool Entity::__puta(const size_t* start, const size_t* size, const unsigned short* values)
{ return _puta_ushort(start, size, values); }

template<> inline
bool Entity::__puta(const size_t* start, const size_t* size, const unsigned int* values)
{ return _puta_uint(start, size, values); }

template<> inline
bool Entity::__puta(const size_t* start, const size_t* size, const long long* values)
{ return _puta_longlong(start, size, values); }

template<> inline
bool Entity::__puta(const size_t* start, const size_t* size, const unsigned long long* values)
{ return _puta_ulonglong(start, size, values); }

template<> inline
bool Entity::__geta(const size_t* start, const size_t* size, signed char* values)
{ return _geta_schar(start, size, values); }

template<> inline
bool Entity::__geta(const size_t* start, const size_t* size, unsigned char* values)
{ return _geta_uchar(start, size, values); }

template<> inline
bool Entity::__geta(const size_t* start, const size_t* size, short* values)
{ return _geta_short(start, size, values); }

template<> inline
bool Entity::__geta(const size_t* start, const size_t* size, int* values)
{ return _geta_int(start, size, values); }

template<> inline
bool Entity::__geta(const size_t* start, const size_t* size, long* values)
{ return _geta_long(start, size, values); }

template<> inline
bool Entity::__geta(const size_t* start, const size_t* size, float* values)
{ return _geta_float(start, size, values); }

template<> inline
bool Entity::__geta(const size_t* start, const size_t* size, double* values)
{ return _geta_double(start, size, values); }

template<> inline
bool Entity::__geta(const size_t* start, const size_t* size, unsigned short* values)
{ return _geta_ushort(start, size, values); }

template<> inline
bool Entity::__geta(const size_t* start, const size_t* size, unsigned int* values)
{ return _geta_uint(start, size, values); }

template<> inline
bool Entity::__geta(const size_t* start, const size_t* size, long long* values)
{ return _geta_longlong(start, size, values); }

template<> inline
bool Entity::__geta(const size_t* start, const size_t* size, unsigned long long* values)
{ return _geta_ulonglong(start, size, values); }

}

#endif // PUML_ENTITY_H
