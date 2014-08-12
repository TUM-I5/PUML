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

#ifndef SIM_MOD_SUITE_H
#define SIM_MOD_SUITE_H

#include <mpi.h>

#include <algorithm>
#include <cassert>
#include <cstring>

#include <SimParasolidKrnl.h>
#include <MeshSim.h>
#include <SimPartitionedMesh.h>
#include <SimError.h>
#include <SimErrorCodes.h>
#include <SimMeshingErrorCodes.h>

#include "utils/logger.h"
#include "utils/progress.h"

#include "MeshGenerator.h"
#include "SimModelerUtil.h"

//forward declare
pAManager SModel_attManager(pModel model);

/**
 * @todo Currently it is not supported to create more than one instance
 *  of this class
 * @todo Maybe add MS_setMaxEntities to limit the number of elements
 */
class SimModSuite : public MeshGenerator
{
private:
	pGModel m_model;

	pAManager m_attMngr;

	pACase m_meshCase;
	pACase m_analysisCase;

	pParMesh m_mesh;
	pMesh m_localMesh;

	/** True of the mesh is already generated */
	bool m_isGenerated;

public:
	SimModSuite(const char* licenseFile = 0L)
		: m_model(0L), m_isGenerated(false)
	{
		_init(licenseFile);
	}

	SimModSuite(const char* cadFile, const char* modFile = 0L,
			const char* licenseFile = 0L,
			const char* meshCaseName = "mesh",
			const char* analysisCaseName = "analysis",
			int enforceSize = 0)
		: m_model(0L), m_isGenerated(false)
	{

		_init(licenseFile);

		generate(cadFile, modFile, meshCaseName, analysisCaseName, enforceSize);
	}

	virtual ~SimModSuite()
	{
		if (m_model) {
			// TODO delete mesh, model, etc
		}

		SimParasolid_stop(1);

		MS_exit();
		Sim_unregisterAllKeys();

		SimPartitionedMesh_stop();
	}

	void generate(const char* cadFile, const char* modFile,
			const char* meshCaseName = "mesh", const char* analysisCaseName = "analysis",
			int forceSize = 0)
	{
		logInfo(PMU_rank()) << "Loading model";
		pNativeModel nativeModel = ParasolidNM_createFromFile(cadFile, 0);

		std::string modelFile;
		if (modFile == 0L) {
			modelFile = cadFile;
			utils::StringUtils::replace(modelFile, "_nat.x_t", ".smd");
		} else
			modelFile = modFile;

		m_model = GM_load(modelFile.c_str(), nativeModel, 0L);
		NM_release(nativeModel);

        // check for model errors
        pPList modelErrors = PList_new();
        if (!GM_isValid(m_model, 0, modelErrors))
                // TODO print more detail about errors
                logError() << "Input model is not valid";
        PList_delete(modelErrors);

		m_attMngr = SModel_attManager(m_model);

		logInfo(PMU_rank()) << "Extracting cases";

		MeshingOptions meshingOptions;
		m_meshCase = MS_newMeshCase(m_model);
		MS_setupSimModelerMeshCase(extractCase(meshCaseName), m_meshCase, &meshingOptions);

		m_analysisCase = extractCase(analysisCaseName);
		pPList children = AttNode_children(m_analysisCase);
		void* iter = 0L;
		while (pANode child = static_cast<pANode>(PList_next(children, &iter)))
			AttCase_setModel(static_cast<pACase>(child), m_model);
		PList_delete(children);

        // create the mesh
        m_mesh = PM_new(0, m_model, PMU_size());

        pProgress prog = Progress_new();
        Progress_setCallback(prog, progressHandler);

        logInfo(PMU_rank()) << "Starting the surface mesher";
        pSurfaceMesher surfaceMesher = SurfaceMesher_new(m_meshCase, m_mesh);
        progressBar.setTotal(26);
        SurfaceMesher_execute(surfaceMesher, prog);

        logInfo(PMU_rank()) << "Starting the volume mesher";
        pVolumeMesher volumeMesher = VolumeMesher_new(m_meshCase, m_mesh);
        VolumeMesher_setEnforceSize(volumeMesher, forceSize);
        progressBar.setTotal(6);
        VolumeMesher_execute(volumeMesher, prog);

        // Cleanup mesher
        SurfaceMesher_delete(surfaceMesher);
        VolumeMesher_delete(volumeMesher);

        Progress_delete(prog);

        // Get basic mesh information

        // Double check if PM_setEntityIds gives the correct ordering
        // Compile without NDEBUG and getElements() and getVertices() will through errors
        int size[4];
        PM_setEntityIds(m_mesh, 1, 1, sthreadDefault, size);
        setNVertices(size[0]);
        PM_setEntityIds(m_mesh, 8, 0, sthreadDefault, size);
        setNElements(size[3]);

        logInfo(PMU_rank()) << "Mesh with" << nElements() << "elements generated";

        // Count slave vertices
        unsigned int numSlaveVertices = 0;
        for (int i = 0; i < PMU_size(); i++) {
        	if (i == PMU_rank())
        		continue;

        	numSlaveVertices += PM_numBdryEntSlave(m_mesh, 0, 0, i);
        }

        m_localMesh = PM_mesh(m_mesh, 0);

        setNLocalVertices(M_numVertices(m_localMesh) - numSlaveVertices);
        setNLocalElements(M_numRegions(m_localMesh));

        init();

        m_isGenerated = true;
	}

	bool isGenerated() const
	{
		return m_isGenerated;
	}

	void getVertices(double* vertices)
	{
        VIter viter = M_vertexIter(m_localMesh);
        while (pVertex v = VIter_next(viter)) {
        	if (!EN_isOwnerProc(v))
        		continue;

        	int localId = EN_id(v) - vertStart(PMU_rank());
        	assert(localId >= 0 && static_cast<unsigned int>(localId) < nLocalVertices());

        	V_coord(v, &vertices[localId*3]);
        }
        VIter_delete(viter);
	}

	void getElements(unsigned int* elements)
	{
		RIter riter = M_regionIter(m_localMesh);
		for (unsigned int i = 0; i < nLocalElements(); i++) {
			pRegion r = RIter_next(riter);
			assert(r);

			int localId = EN_id(r) - elemStart(PMU_rank());
			assert(localId >= 0 && static_cast<unsigned int>(localId) < nLocalElements());

			pPList vertices = R_vertices(r, 1); // 1 should be the correct ordering ...
			// TODO assuming tetrahedra
			for (unsigned int j = 0; j < 4; j++)
				elements[localId*4+j] = EN_id(static_cast<pVertex>(PList_item(vertices, j)));
			PList_delete(vertices);
		}
		RIter_delete(riter);
	}

	void getGroups(unsigned int* groups)
	{
		for (unsigned int i = 0; i < nLocalElements(); i++)
			// TODO get the group from the model region
			groups[i] = 0;
	}

	void getBoundaries(unsigned int* boundaries)
	{
		AttCase_associate(m_analysisCase, 0L);

		FIter fiter = M_faceIter(m_localMesh);
		while (pFace face = FIter_next(fiter)) {
			if (F_whatInType(face) != Gface)
				continue;

			pGFace modelFace = static_cast<pGFace>(F_whatIn(face));
			pAttribute attr = GEN_attrib(modelFace, "boundaryCondition");
			if (attr) {
				char* image = Attribute_imageClass(attr);
				unsigned int boundary = parseBoundary(image);
				Sim_deleteString(image);

				for (int i = 0; i <= 1; i++) {
					// Set boundary for both sides

					pRegion r = F_region(face, i);
					if (!r)
						continue;

					for (unsigned int j = 0; j < 4; j++) {
						pFace f = R_face(r, j);

						if (f == face) {
							int localId = EN_id(r) - elemStart(PMU_rank());
							assert(localId >= 0
									&& static_cast<unsigned int>(localId) < nLocalElements());

							boundaries[localId*4 + FACE2INTERNAL[j]] = boundary;

							break;
						}
					}
				}
			}
		}
		FIter_delete(fiter);

		AttCase_unassociate(m_analysisCase);
	}

	/**
	 * Provides direct access to the SimModSuite model. Should be used with caution.
	 */
	pGModel model()
	{
		return m_model;
	}

	/**
	 * Provides direct access to the SimModSuite mesh. Should be used with caution.
	 */
	pParMesh mesh()
	{
		return m_mesh;
	}

	/**
	 * Provides direct access to the SimModSuite analysis case. Should be used with caution.
	 */
	pACase analysisCase()
	{
		return m_analysisCase;
	}

private:
	pACase extractCase(const char* name)
	{
		pACase acase = AMAN_findCase(m_attMngr, name);
		if (!acase)
			logError() << "Case" << std::string(name) << "not found.";

		AttCase_setModel(acase, m_model);

		return acase;
	}

public:
	/**
	 * @todo Make this private as soon as APF does no longer need it
	 */
	static unsigned int parseBoundary(const char* boundaryCondition)
	{
		if (strcmp(boundaryCondition, "freeSurface") == 0)
			return 1;
		if (strcmp(boundaryCondition, "dynamicRupture") == 0)
			return 3;
		if (strcmp(boundaryCondition, "absorbing") == 0)
			return 5;

		logError() << "Unknown boundary condition" << boundaryCondition;
		return -1;
	}

private:
	static void _init(const char* licenseFile)
	{
		SimPartitionedMesh_start(0L, 0L);
		Sim_readLicenseFile(licenseFile);
		MS_init();
		SimParasolid_start(1);

		Sim_setMessageHandler(messageHandler);
	}

	static void messageHandler(int type, const char* msg)
	{
		switch (type) {
		case Sim_InfoMsg:
			// Show sim info messages as debug messages
			logDebug(PMU_rank()) << "SimModeler:" << msg;
			break;
		case Sim_DebugMsg:
			// Ignore sim debug messages
			break;
		case Sim_WarningMsg:
			logWarning(PMU_rank()) << "SimModeler:" << msg;
			break;
		case Sim_ErrorMsg:
			// Use warning because error will abort the program
			logWarning() << "SimModeler:" << msg;
			break;
		}
	}

	static void progressHandler(const char* what, int level, int startVal, int endVal, int currentVal, void *ignore)
	{
		if (PMU_rank() != 0 || level >= 2)
			return;

		switch (level) {
		case 0:
			if (currentVal == -2)
				progressBar.update(0);
			else
				progressBar.clear();
			break;
		case 1:
			if (currentVal == 0)
				progressBar.update();
			else
				progressBar.increment();
			break;
		}

		logDebug() << what << level << startVal << endVal << currentVal;
	}

public:
	static const unsigned int FACE2INTERNAL[];

private:
	static utils::Progress progressBar;

};

#endif // SIM_MOD_SUITE_H
