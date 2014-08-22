/**
 * @file
 *  This file is part of PUML
 *
 *  For conditions of distribution and use, please see the copyright
 *  notice in the file 'COPYING' at the root directory of this package
 *  and the copyright notice at https://github.com/TUM-I5/PUML
 *
 * @copyright 2014 Technische Universitaet Muenchen
 * @author Sebastian Rettenberger <rettenbs@in.tum.de>
 */

#ifndef INPUT_APF_H
#define INPUT_APF_H

#include <limits>

#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <gmi.h>
#include <maMesh.h>
#include <PCU.h>

#include "utils/logger.h"

class Apf
{
private:
	/** Number of processes */
	int m_nProcs;

	/** APF mesh */
	apf::Mesh2* m_mesh;

	/** Global numbering for the vertices */
	apf::GlobalNumbering* m_vertexNum;

	/** Global number of the elements */
	apf::GlobalNumbering* m_elementNum;

	/** Tag for the boundary condition */
	apf::MeshTag* m_boundaryTag;

	/** Tag for the group */
	apf::MeshTag* m_groupTag;

	unsigned int m_nVertices;
	unsigned int m_nElements;

	/** Number of local vertices */
	unsigned int m_nLocalVertices;
	/** Number of local elements */
	unsigned int m_nLocalElements;

	/** The global id of the first vertex on each process */
	unsigned int* m_vertStart;
	/** The global id of the first element on each process */
	unsigned int* m_elemStart;

public:
	/**
	 * @param mesh An apf mesh. The mesh will be automatically
	 *  destoryed when this instance ins deleted.
	 */
	Apf(apf::Mesh2* mesh)
		: m_mesh(mesh)
	{
		int rank = PCU_Comm_Self();
		m_nProcs = PCU_Comm_Peers();

		// Check mesh
		if (alignMdsMatches(m_mesh))
			logWarning() << "Fixed misaligned matches";
		mesh->verify();

		// Create numberings
		m_vertexNum = apf::makeGlobal(
				apf::numberOwnedDimension(mesh, "vertices", 0));
		apf::synchronize(m_vertexNum);
		m_elementNum = apf::makeGlobal(
				apf::numberOwnedDimension(mesh, "elements", 3));

		// Get tags
		m_boundaryTag = mesh->findTag("boundary condition");
		m_groupTag = mesh->findTag("group");

		// Get local size
		m_nLocalVertices = apf::countOwned(mesh, 0);
		m_nLocalElements = apf::countOwned(mesh, 3);

		// Get the total number of vertices/elements
		unsigned int size[2] = {m_nLocalVertices, m_nLocalElements};
		MPI_Allreduce(MPI_IN_PLACE, size, 2, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

		m_nVertices = size[0];
		m_nElements = size[1];

		// Initialize mapping arrays
        unsigned int start[2] = {m_nLocalVertices, m_nLocalElements};
        MPI_Scan(MPI_IN_PLACE, start, 2, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

        m_vertStart = new unsigned int[m_nProcs];
        m_vertStart[rank] = start[0] - m_nLocalVertices;
        MPI_Allgather(MPI_IN_PLACE, 1, MPI_UNSIGNED, m_vertStart, 1, MPI_UNSIGNED,
        		MPI_COMM_WORLD);

        m_elemStart = new unsigned int[m_nProcs];
        m_elemStart[rank] = start[1] - m_nLocalElements;
        MPI_Allgather(MPI_IN_PLACE, 1, MPI_UNSIGNED, m_elemStart, 1, MPI_UNSIGNED,
        		MPI_COMM_WORLD);

		// Compute min insphere radius
		double min = std::numeric_limits<double>::max();
		apf::MeshIterator* eiter = m_mesh->begin(3);
		while (apf::MeshEntity* element = m_mesh->iterate(eiter)) {
			min = std::min(min, ma::getInsphere(m_mesh, element));
		}
		MPI_Reduce((PCU_Comm_Self() == 0 ? MPI_IN_PLACE : &min),
				&min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
		logInfo(PCU_Comm_Self()) << "Minimum insphere found:" << min;
	}

	virtual ~Apf()
	{
		apf::destroyGlobalNumbering(m_vertexNum);
		apf::destroyGlobalNumbering(m_elementNum);

		if (m_boundaryTag)
			m_mesh->destroyTag(m_boundaryTag);
		if (m_groupTag)
			m_mesh->destroyTag(m_groupTag);

		m_mesh->destroyNative();
		apf::destroyMesh(m_mesh);

		delete [] m_vertStart;
		delete [] m_elemStart;
	}

	unsigned int nVertices() const
	{
		return m_nVertices;
	}

	unsigned int nElements() const
	{
		return m_nElements;
	}

	/**
	 * Number of vertices this process is responsible for
	 */
	unsigned int nLocalVertices() const
	{
		return m_nLocalVertices;
	}

	unsigned int nLocalElements() const
	{
		return m_nLocalElements;
	}

	unsigned int vertStart(int rank) const
	{
		return m_vertStart[rank];
	}

	int rankOfVert(unsigned int vertex) const
	{
		return std::upper_bound(m_vertStart, m_vertStart+m_nProcs, vertex) - m_vertStart - 1;
	}

	unsigned int posOfVert(unsigned int vertex) const
	{
		return vertex - m_vertStart[rankOfVert(vertex)];
	}

	unsigned int elemStart(int rank) const
	{
		return m_elemStart[rank];
	}

	int rankOfElem(unsigned int element) const
	{
		return std::upper_bound(m_elemStart, m_elemStart+m_nProcs, element) - m_elemStart - 1;
	}

	unsigned int posOfElem(unsigned int element) const
	{
		return element - m_elemStart[rankOfElem(element)];
	}

	void getVertices(double* vertices)
	{
		apf::MeshIterator* viter = m_mesh->begin(0);
		while (apf::MeshEntity* vertex = m_mesh->iterate(viter)) {
			if (!m_mesh->isOwned(vertex))
				continue;

			long localId = apf::getNumber(m_vertexNum, apf::Node(vertex, 0))
					- vertStart(PCU_Comm_Self());
			assert(localId >= 0 && static_cast<unsigned int>(localId) < nLocalVertices());

			apf::Vector3 point;
			m_mesh->getPoint(vertex, 0, point);
			point.toArray(&vertices[localId*3]);
		}
		m_mesh->end(viter);
	}

	void getElements(unsigned int* elements)
	{
		apf::MeshIterator* eiter = m_mesh->begin(3);
		while (apf::MeshEntity* element = m_mesh->iterate(eiter)) {
			assert(m_mesh->isOwned(element));

			long localId = apf::getNumber(m_elementNum, apf::Node(element, 0))
					- elemStart(PCU_Comm_Self());
			assert(localId >= 0 && static_cast<unsigned int>(localId) < nLocalElements());

			apf::Adjacent adjacent;
			m_mesh->getAdjacent(element, 0, adjacent);
			assert(adjacent.getSize() == 4);
			// TODO check correct ordering of vertices

			for (size_t i = 0; i < 4; i++)
				elements[localId*4+i] = apf::getNumber(m_vertexNum, apf::Node(adjacent[i], 0));
		}
		m_mesh->end(eiter);
	}

	void getGroups(unsigned int* groups)
	{
		if (!m_groupTag) {
			memset(groups, 0, nLocalElements()*sizeof(unsigned int));
			return;
		}

		apf::MeshIterator* it = m_mesh->begin(3);
		while (apf::MeshEntity* element = m_mesh->iterate(it)) {
			assert(m_mesh->isOwned(element));
			assert(m_mesh->hasTag(element, m_groupTag));

			long localId = apf::getNumber(m_elementNum, apf::Node(element, 0))
					- elemStart(PCU_Comm_Self());
			assert(localId >= 0 && static_cast<unsigned int>(localId) < nLocalElements());

			int group;
			m_mesh->getIntTag(element, m_groupTag, &group);
			groups[localId] = group;
		}
		m_mesh->end(it);
	}

	void getBoundaries(unsigned int* boundaries)
	{
		if (!m_boundaryTag)
			return;

		apf::MeshIterator* it = m_mesh->begin(3);
		while (apf::MeshEntity* element = m_mesh->iterate(it)) {
			assert(m_mesh->isOwned(element));

			long localId = apf::getNumber(m_elementNum, apf::Node(element, 0))
					- elemStart(PCU_Comm_Self());
			assert(localId >= 0 && static_cast<unsigned int>(localId) < nLocalElements());

			apf::Adjacent adjacent;
			m_mesh->getAdjacent(element, 2, adjacent);
			assert(adjacent.getSize() == 4);

			for (size_t i = 0; i < 4; i++) {
				if (!m_mesh->hasTag(adjacent[i], m_boundaryTag))
					continue;

				int tag;
				m_mesh->getIntTag(adjacent[i], m_boundaryTag, &tag);
				boundaries[localId*4 + FACE2INTERNAL[i]] = tag;
			}
		}
		m_mesh->end(it);
	}

	void write(const char* mesh, const char* model = 0L)
	{
		logInfo(PCU_Comm_Self()) << "Writing native APF mesh";
		m_mesh->writeNative(mesh);
		if (model)
			gmi_write_dmg(m_mesh->getModel(), model);
	}

private:
	static const unsigned int FACE2INTERNAL[];
};

#endif // INPUT_APF_H
