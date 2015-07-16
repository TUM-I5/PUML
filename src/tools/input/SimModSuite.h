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
 * @author Johannes Klicpera <johannes.klicpera@tum.de>
 */

#ifndef SIM_MOD_SUITE_H
#define SIM_MOD_SUITE_H

#include <MeshSim.h>
#include <SimPartitionedMesh.h>

#include "utils/progress.h"

#include "MeshInput.h"

//forward declare
pAManager SModel_attManager(pModel model);

/**
 * @todo Currently it is not supported to create more than one instance
 *  of this class
 * @todo Maybe add MS_setMaxEntities to limit the number of elements
 */
class SimModSuite : public MeshInput
{
private:
	pGModel m_model;
	pParMesh m_simMesh;

	/** Enable Simmetrix logging file */
	bool m_log;

public:
	SimModSuite(const char* modFile, const char* cadFile = 0L,
			const char* licenseFile = 0L,
			const char* meshCaseName = "mesh",
			const char* analysisCaseName = "analysis",
			int enforceSize = 0,
			const char* stl_ParFile = 0L,
			const bool probe_faces = false,
			const char* logFile = 0L);
	virtual ~SimModSuite();

private:
	pACase extractCase(pAManager attMngr, const char* name);
	static unsigned int parseBoundary(const char* boundaryCondition);
	static void messageHandler(int type, const char* msg);
	static void progressHandler(const char* what, int level, int startVal, int endVal, int currentVal, void *ignore);
	static utils::Progress progressBar;

	void loadCAD(const char* modFile, const char* cadFile);
	void loadSTL(const char *filename);
	void extractCases(pGModel m_model, pACase &meshCase, const char *meshCaseName, pACase &analysisCase, const char *analysisCaseName);
	void setCases(pGModel model, pACase &meshCase, pACase &analysisCase, const char* stl_ParFile);
	void analyse_mesh();
	void probeFaceCoords(pGModel model);
};

#endif // SIM_MOD_SUITE_H
