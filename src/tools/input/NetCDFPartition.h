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

#include <algorithm>
#include <cassert>
#include <cstring>
#include <map>
#include <vector>

#include "utils/logger.h"

/**
 * Describes one partition (required for reading netCDF meshes)
 *
 * @todo Global/local names are probably wrong
 */
class Partition
{
private:
	struct Face {
		unsigned int element;
		unsigned int face;
	};

private:
	int m_rank;

	unsigned int m_nElements;
	unsigned int m_nVertices;
	unsigned int m_nLocalVertices;

	int* m_elements;
	int* m_neighborRanks;
	int* m_neighborOrientation;
	int* m_mpiIndices;
	int* m_boundaries;
	int* m_groups;

	double* m_vertices;

	int* m_globalElements;
	double* m_globalVertices;

	int* m_vertexOffset;
	/** Vertices where this partition is responsible */
	bool* m_isLocalVertices;

	std::map<int, std::vector<Face>> m_neighborLists;

	/** Vertex map */
	std::map<int, int> m_local2global;

public:
	Partition()
		: m_rank(0),
		  m_nElements(0), m_nVertices(0), m_nLocalVertices(0),
		  m_elements(0L), m_neighborRanks(0L),
		  m_neighborOrientation(0L), m_mpiIndices(0L),
		  m_boundaries(0L), m_groups(0L),
		  m_vertices(0L),
		  m_globalElements(0L), m_globalVertices(0L),
		  m_vertexOffset(0L), m_isLocalVertices(0L)
	{}

	~Partition()
	{
		delete [] m_elements;
		delete [] m_neighborRanks;
		delete [] m_neighborOrientation;
		delete [] m_mpiIndices;
		delete [] m_boundaries;
		delete [] m_groups;

		delete [] m_vertices;

		delete [] m_globalElements;
		delete [] m_globalVertices;

		delete [] m_vertexOffset;
		delete [] m_isLocalVertices;
	}

	void setRank(int rank)
	{
		m_rank = rank;
	}

	void setElemSize(unsigned int nElements)
	{
		if (m_nElements != 0)
			return;

		m_nElements = nElements;

		m_elements = new int[nElements*4];
		m_neighborRanks = new int[nElements*4];
		m_neighborOrientation = new int[nElements*4];
		m_mpiIndices = new int[nElements*4];
		m_boundaries = new int[nElements*4];
		m_groups = new int[nElements];

		m_globalElements = new int[nElements*4];
		for (unsigned int i = 0; i < nElements*4; i++)
			m_globalElements[i] = -1;
	}

	void setVrtxSize(unsigned int nVertices)
	{
		if (m_nVertices != 0)
			return;

		m_nVertices = nVertices;

		m_vertices = new double[nVertices*3];
		m_vertexOffset = new int[nVertices];
		m_isLocalVertices = new bool[nVertices];
		std::fill_n(m_isLocalVertices, nVertices, true);
	}

	void computeLocalVertices()
	{
		for (unsigned int i = 0; i < m_nElements*4; i++) {
			if (m_neighborRanks[i] < m_rank) {
				unsigned int face = i % 4;
				unsigned int element = i / 4;

				for (unsigned int j = 0; j < 3; j++)
					m_isLocalVertices[m_elements[element*4 + FACE2VERTICES[face][j]]] = false;
			}
		}

		m_vertexOffset[0] = 0;
		for (unsigned int i = 1; i < m_nVertices; i++) {
			if (m_isLocalVertices[i-1]) {
				m_vertexOffset[i] = m_vertexOffset[i-1] + 1;
				m_nLocalVertices++;
			} else
				m_vertexOffset[i] = m_vertexOffset[i-1];
		}
		if (m_isLocalVertices[m_nVertices-1])
			m_nLocalVertices++;
	}

	void convertLocalVertices(unsigned int startId)
	{
		for (unsigned int i = 0; i < m_nElements*4; i++)
			if (m_isLocalVertices[m_elements[i]])
				m_globalElements[i] = m_vertexOffset[m_elements[i]] + startId;
	}

	void buildVertexRecvLists()
	{
		for (unsigned int i = 0; i < m_nElements*4; i++) {
			if (m_neighborRanks[i] < m_rank) {
				// Compute neighbor list for later usage
				unsigned int face = i % 4;
				unsigned int element = i / 4;

				std::vector<Face> &list = m_neighborLists[m_neighborRanks[i]];
				if (list.size() <= static_cast<size_t>(m_mpiIndices[i]))
					list.resize(m_mpiIndices[i]+1);
				list[m_mpiIndices[i]].element = element;
				list[m_mpiIndices[i]].face = face;
			}
		}
	}

	void buildVertexSendLists(std::map<int, std::vector<int>> &globalVertexIds)
	{
		for (unsigned int i = 0; i < m_nElements*4; i++) {
			if (m_neighborRanks[i] > m_rank) {
				unsigned int face = i % 4;
				unsigned int element = i / 4;

				std::vector<int> &list = globalVertexIds[m_neighborRanks[i]];
				if (list.size() <= static_cast<size_t>(m_mpiIndices[i])*3)
					list.resize(m_mpiIndices[i]*3+3);
				for (unsigned int j = 0; j < 3; j++) {
					list[m_mpiIndices[i]*3+j] = m_globalElements[element*4 + FACE2VERTICES[face][j]];
				}
			}
		}
	}

	bool buildLocal2GlobalMap(int remoteRank, int* vertexIds)
	{
		std::vector<Face> &list = m_neighborLists[remoteRank];
		assert(list.size() > 0);

		bool done = true;

		for (unsigned int i = 0; i < list.size(); i++) {
			int orientation = m_neighborOrientation[list[i].element*4 + list[i].face];

			int localVertices[3] = {
					m_elements[list[i].element*4 + FACE2VERTICES[list[i].face][0]],
					m_elements[list[i].element*4 + FACE2VERTICES[list[i].face][1]],
					m_elements[list[i].element*4 + FACE2VERTICES[list[i].face][2]]
			};

			int globalVertices[3] = {
					vertexIds[i*3 + orientation],
					vertexIds[i*3 + (orientation+2) % 3],
					vertexIds[i*3 + (orientation+1) % 3]
			};

			for (unsigned int j = 0; j < 3; j++) {
				if (globalVertices[j] < 0) {
					// We got a wrong value, update the elements and resend them
					done = false;
					continue;
				}

				std::map<int, int>::const_iterator it = m_local2global.find(localVertices[j]);

				if (it != m_local2global.end()) {
					if (it->second != globalVertices[j])
						logError() << "Invalid mapping detected";
				} else
					m_local2global[localVertices[j]] = globalVertices[j];
			}
		}

		return done;
	}

	void applyLocal2GlobalMap()
	{
		for (unsigned int i = 0; i < m_nElements*4; i++) {
			std::map<int, int>::const_iterator it = m_local2global.find(m_elements[i]);

			if (it != m_local2global.end())
				m_globalElements[i] = it->second;
		}
	}

	void extractGlobalVertices()
	{
		m_globalVertices = new double[m_nLocalVertices*3];

		unsigned int pos = 0;
		for (unsigned int i = 0; i < m_nVertices; i++) {
			if (m_isLocalVertices[i]) {
				memcpy(&m_globalVertices[pos*3], &m_vertices[i*3], 3*sizeof(double));
				pos++;
			}
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

	int rank() const
	{
		return m_rank;
	}

	unsigned int nElements() const
	{
		return m_nElements;
	}

	unsigned int nVertices() const
	{
		return m_nVertices;
	}

	unsigned int nLocalVertices() const
	{
		return m_nLocalVertices;
	}

	int* elements()
	{
		return m_elements;
	}

	int* neighborRanks()
	{
		return m_neighborRanks;
	}

	int* neighborOrientation()
	{
		return m_neighborOrientation;
	}

	int* mpiIndices()
	{
		return m_mpiIndices;
	}

	int* boundaries()
	{
		return m_boundaries;
	}

	int* groups()
	{
		return m_groups;
	}

	double* vertices()
	{
		return m_vertices;
	}

	const int* globalElements() const
	{
		return m_globalElements;
	}

	const double* globalVertices() const
	{
		return m_globalVertices;
	}

private:
	const static int FACE2VERTICES[4][3];
	const static int INTERNAL2EX_ORDER[4];
};

#endif // PARTITION_H
