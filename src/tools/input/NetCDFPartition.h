/**
 * @file
 *  This file is part of PUML
 *
 *  For conditions of distribution and use, please see the copyright
 *  notice in the file 'COPYING' at the root directory of this package
 *  and the copyright notice at https://github.com/TUM-I5/PUML
 *
 * @copyright 2015 Technische Universitaet Muenchen
 * @author Sebastian Rettenberger <rettenbs@in.tum.de>
 */

#ifndef PARTITION_H
#define PARTITION_H

#include <cstring>

/**
 * Describes one partition (required for reading netCDF meshes)
 */
class Partition
{
private:
	unsigned int m_nElements;
	unsigned int m_nVertices;

	int* m_elements;
	double* m_vertices;
	int* m_boundaries;
	int* m_groups;

public:
	Partition()
		: m_nElements(0), m_nVertices(0),
		  m_elements(0L), m_vertices(0L),
		  m_boundaries(0L), m_groups(0L)
	{}

	~Partition()
	{
		delete [] m_elements;
		delete [] m_vertices;
		delete [] m_boundaries;
		delete [] m_groups;
	}

	void setElemSize(unsigned int nElements)
	{
		if (m_nElements != 0)
			return;

		m_nElements = nElements;

		m_elements = new int[nElements*4];
		m_boundaries = new int[nElements*4];
		m_groups = new int[nElements];
	}

	void setVrtxSize(unsigned int nVertices)
	{
		if (m_nVertices != 0)
			return;

		m_nVertices = nVertices;

		m_vertices = new double[nVertices*3];
	}
	
	void setEqualGroups(int group){	
	  for(unsigned int i = 0 ; i < m_nElements; i++){
	    m_groups[i] = group;
	  }
	}

	void convertBoundary()
	{
		int ncBoundaries[4];

		for (unsigned int i = 0; i < m_nElements*4; i += 4) {
			memcpy(ncBoundaries, &m_boundaries[i], 4*sizeof(int));
			for (unsigned int j = 0; j < 4; j++)
				m_boundaries[i+j] = ncBoundaries[INTERNAL2EX_ORDER[j]];
		}
	}

	unsigned int nElements() const
	{
		return m_nElements;
	}

	unsigned int nVertices() const
	{
		return m_nVertices;
	}

	int* elements()
	{
		return m_elements;
	}

	double* vertices()
	{
		return m_vertices;
	}

	int* boundaries()
	{
		return m_boundaries;
	}

	int* groups()
	{
		return m_groups;
	}

private:
	const static int INTERNAL2EX_ORDER[4];
};

#endif // PARTITION_H
