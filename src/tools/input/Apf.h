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

#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <apfSIM.h>
#include <apfShape.h>
#include <gmi_sim.h>
#include <PCU.h>

#include "utils/logger.h"

#include "MeshGenerator.h"
#include "SimModSuite.h"

class Apf : public MeshGenerator
{
private:
	/** APF mesh */
	apf::Mesh2* m_mesh;

	/** Global numbering for the vertices */
	apf::GlobalNumbering* m_vertexNum;

	/** Global number of the elements */
	apf::GlobalNumbering* m_elementNum;

	/** Tag for the boundary condition */
	apf::MeshTag* m_boundaryTag;

public:
	Apf()
		: m_mesh(0L)
	{
		_init();
	}

	Apf(SimModSuite &meshSource)
		: m_mesh(0L)
	{
		_init();
		load(meshSource);
	}

	virtual ~Apf()
	{
		if (m_mesh) {
			apf::destroyGlobalNumbering(m_vertexNum);
			apf::destroyGlobalNumbering(m_elementNum);

			m_mesh->destroyTag(m_boundaryTag);

			m_mesh->destroyNative();
			apf::destroyMesh(m_mesh);
		}

		PCU_Comm_Free();
	}

	/**
	 * Initializes APF data structures. Has to be called after the mesh is
	 * generated.
	 */
	void load(SimModSuite &meshSource)
	{
		if (!meshSource.isGenerated())
			logError() << "Can not initialize APF, mesh not generated.";

		apf::Mesh* mesh = apf::createMesh(meshSource.mesh());
		gmi_register_sim();
		gmi_model* model = gmi_import_sim(meshSource.model());

		logInfo(PCU_Comm_Self()) << "Converting data structure";
		m_mesh = apf::createMdsMesh(model, mesh);
		apf::destroyMesh(mesh);

		if (alignMdsMatches(m_mesh))
			logWarning() << "Fixed misaligned matches";
		m_mesh->verify();

		m_vertexNum = apf::makeGlobal(
				apf::numberOwnedDimension(m_mesh, "vertices", 0));
		apf::synchronize(m_vertexNum);
		m_elementNum = apf::makeGlobal(
				apf::numberOwnedDimension(m_mesh, "elements", 3));

		setNLocalVertices(apf::countOwned(m_mesh, 0));
		setNLocalElements(apf::countOwned(m_mesh, 3));

		// Get the total number of vertices/elements
		unsigned int size[2] = {nLocalVertices(), nLocalElements()};
		MPI_Allreduce(MPI_IN_PLACE, size, 2, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

		setNVertices(size[0]);
		setNElements(size[1]);

		init();

		// Set the boundary conditions from the geometric model
		AttCase_associate(meshSource.analysisCase(), 0L);
		m_boundaryTag = m_mesh->createIntTag("boundary condition", 1);
		apf::MeshIterator* fiter = m_mesh->begin(2);
		while (apf::MeshEntity* face = m_mesh->iterate(fiter)) {
			apf::ModelEntity* modelFace = m_mesh->toModel(face);
			if (m_mesh->getModelType(modelFace) != 2)
				continue;

			pGEntity simFace = reinterpret_cast<pGEntity>(modelFace);

			pAttribute attr = GEN_attrib(simFace, "boundaryCondition");
			if (attr) {
				char* image = Attribute_imageClass(attr);
				int boundary = SimModSuite::parseBoundary(image);
				Sim_deleteString(image);

				m_mesh->setIntTag(face, m_boundaryTag, &boundary);
			}

		}
		m_mesh->end(fiter);
		AttCase_unassociate(meshSource.analysisCase());
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
		for (unsigned int i = 0; i < nLocalElements(); i++)
			// TODO get the group from the model region
			groups[i] = 0;
	}

	void getBoundaries(unsigned int* boundaries)
	{
		apf::MeshIterator* eiter = m_mesh->begin(3);
		while (apf::MeshEntity* element = m_mesh->iterate(eiter)) {
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
				// TODO do not use simmodsuite here
				m_mesh->getIntTag(adjacent[i], m_boundaryTag, &tag);
				boundaries[localId*4+SimModSuite::FACE2INTERNAL[i]] = tag;
			}
		}
		m_mesh->end(eiter);
	}

private:
	/**
	 * Common initialization for all constructors
	 */
	static void _init()
	{
		PCU_Comm_Init();
		PCU_Protect();
	}
};

#endif // INPUT_APF_H
