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

 #include "SimModSuite.h"

 #include <mpi.h>

 #include <algorithm>
 #include <cassert>
 #include <cstring>

 #include <apf.h>
 #include <apfMDS.h>
 #include <apfMesh2.h>
 #include <apfSIM.h>
 #include <gmi_sim.h>

#include <SimDiscrete.h>
#include <SimParasolidKrnl.h>
#include <SimError.h>
#include <SimErrorCodes.h>
#include <SimMeshingErrorCodes.h>

#include "utils/logger.h"
#include "utils/path.h"

#include "SimModelerUtil.h"

//forward declare
pAManager SModel_attManager(pModel model);

utils::Progress SimModSuite::progressBar;

SimModSuite::SimModSuite(const char* modFile, const char* cadFile,
        const char* licenseFile,
        const char* meshCaseName,
        const char* analysisCaseName,
        int enforceSize,
        const bool stlInput,
        const char* logFile)
{

    // Init SimModSuite
    SimPartitionedMesh_start(0L, 0L);
    if (logFile) {
        m_log = true;
        Sim_logOn(logFile);
    } else
        m_log = false;
    Sim_readLicenseFile(licenseFile);
    MS_init();
    SimDiscrete_start(0);
    SimParasolid_start(1);
    Sim_setMessageHandler(messageHandler);

    // Load file
    logInfo(PMU_rank()) << "Loading model";

    if(stlInput) {
        loadSTL(modFile);
    } else {
        loadCAD(modFile, cadFile);
    }

    // Extract cases
    pACase meshCase, analysisCase;
    if(stlInput) {
        setCases(m_model, meshCase, analysisCase);
    } else {
        extractCases(m_model, meshCase, meshCaseName, analysisCase, analysisCaseName);
    }

    // create the mesh
    m_simMesh = PM_new(0, m_model, PMU_size());

    pProgress prog = Progress_new();
    Progress_setCallback(prog, progressHandler);

    logInfo(PMU_rank()) << "Starting the surface mesher";
    pSurfaceMesher surfaceMesher = SurfaceMesher_new(meshCase, m_simMesh);
    progressBar.setTotal(26);
    SurfaceMesher_execute(surfaceMesher, prog);

    logInfo(PMU_rank()) << "Starting the volume mesher";
    pVolumeMesher volumeMesher = VolumeMesher_new(meshCase, m_simMesh);
    VolumeMesher_setEnforceSize(volumeMesher, enforceSize);
    progressBar.setTotal(6);
    VolumeMesher_execute(volumeMesher, prog);

    // Cleanup mesher
    SurfaceMesher_delete(surfaceMesher);
    VolumeMesher_delete(volumeMesher);

    Progress_delete(prog);

    // Convert to APF mesh
    apf::Mesh* tmpMesh = apf::createMesh(m_simMesh);
    gmi_register_sim();
    gmi_model* model = gmi_import_sim(m_model);

    logInfo(PMU_rank()) << "Converting mesh to APF";
    m_mesh = apf::createMdsMesh(model, tmpMesh);
    apf::destroyMesh(tmpMesh);

    // Set the boundary conditions from the geometric model
    logInfo(PMU_rank()) << "Setting boundary conditions";
    AttCase_associate(analysisCase, 0L);
    apf::MeshTag* boundaryTag = m_mesh->createIntTag("boundary condition", 1);
    apf::MeshIterator* it = m_mesh->begin(2);
    while (apf::MeshEntity* face = m_mesh->iterate(it)) {
        apf::ModelEntity* modelFace = m_mesh->toModel(face);
        if (m_mesh->getModelType(modelFace) != 2)
            continue;

        pGEntity simFace = reinterpret_cast<pGEntity>(modelFace);

        pAttribute attr = GEN_attrib(simFace, "boundaryCondition");
        if (attr) {
            char* image = Attribute_imageClass(attr);
            int boundary = parseBoundary(image);
            Sim_deleteString(image);

            m_mesh->setIntTag(face, boundaryTag, &boundary);
        }

    }
    m_mesh->end(it);
    AttCase_unassociate(analysisCase);

    // Delete cases
    MS_deleteMeshCase(meshCase);
    MS_deleteMeshCase(analysisCase);
}

SimModSuite::~SimModSuite()
{
    M_release(m_simMesh);
    // TODO we can delete the model here because it is still
    // connected to the mesh
    //GM_release(m_model);

    // Finalize SimModSuite
    SimParasolid_stop(1);
    SimDiscrete_stop(0);
    MS_exit();
    Sim_unregisterAllKeys();
    if (m_log)
        Sim_logOff();
    SimPartitionedMesh_stop();
}

pACase SimModSuite::extractCase(pAManager attMngr, const char* name)
{
    pACase acase = AMAN_findCase(attMngr, name);
    if (!acase)
        logError() << "Case" << std::string(name) << "not found.";

    AttCase_setModel(acase, m_model);

    return acase;
}

unsigned int SimModSuite::parseBoundary(const char* boundaryCondition)
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

void SimModSuite::messageHandler(int type, const char* msg)
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

void SimModSuite::progressHandler(const char* what, int level, int startVal, int endVal, int currentVal, void *ignore)
{
	if (PMU_rank() != 0)
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

void SimModSuite::extractCases(pGModel m_model, pACase &meshCase, const char *meshCaseName, pACase &analysisCase, const char *analysisCaseName) {
    logInfo(PMU_rank()) << "Extracting cases";
    pAManager attMngr = SModel_attManager(m_model);

    MeshingOptions meshingOptions;
    meshCase = MS_newMeshCase(m_model);
    MS_setupSimModelerMeshCase(extractCase(attMngr, meshCaseName),
            meshCase, &meshingOptions);

    analysisCase = extractCase(attMngr, analysisCaseName);
    pPList children = AttNode_children(analysisCase);
    void* iter = 0L;
    while (pANode child = static_cast<pANode>(PList_next(children, &iter)))
        AttCase_setModel(static_cast<pACase>(child), m_model);
    PList_delete(children);
}

void SimModSuite::setCases(pGModel model, pACase &meshCase, pACase &analysisCase) {

    logInfo(PMU_rank()) << "Setting cases";
    // ------------------------------ Set boundary conditions ------------------------------

    // Create a new manager. The manager is responsible for creating
    // new attributes, saving/retrieving to/from file, and to look up
    // cases, AttModels etc.
    pAManager attMngr = AMAN_new();

    // Create a case. A case serves as a grouping mechanism. It will (should)
    // contain all the attributes that make up a complete analysis of
    // some kind, although the system has no way of verifying this.
    analysisCase = AMAN_newCase(attMngr,"analysis","",(pModel)model);

    // Now we create an attribute information node. This can be seen
    // as an "attribute generator object", as we will later on create an
    // actual attribute for each geometric face from this attribute
    // information node. We name the attribute T1, and give it the
    // information type "boundaryCondition".
    pAttInfoVoid iSurf = AMAN_newAttInfoVoid(attMngr,"BC","boundaryCondition");
    AttNode_setImageClass((pANode)iSurf,"freeSurface");
    pAttInfoVoid iDynRup = AMAN_newAttInfoVoid(attMngr,"BC","boundaryCondition");
    AttNode_setImageClass((pANode)iDynRup,"dynamicRupture");
    pAttInfoVoid iAbsorb = AMAN_newAttInfoVoid(attMngr,"BC","boundaryCondition");
    AttNode_setImageClass((pANode)iAbsorb,"absorbing");

    // We need to add the Attribute Information Nodes to the case
    AttCase_addNode(analysisCase,(pANode)iSurf);
    AttCase_addNode(analysisCase,(pANode)iDynRup);
    AttCase_addNode(analysisCase,(pANode)iAbsorb);

    // We need to associate the Attribute Information Nodes with the
    // model. To do that we have to create a Model Association for the case
    pModelAssoc aSurf = AttCase_newModelAssoc(analysisCase,(pANode)iSurf);
    pModelAssoc aDynRup = AttCase_newModelAssoc(analysisCase,(pANode)iDynRup);
    pModelAssoc aAbsorb = AttCase_newModelAssoc(analysisCase,(pANode)iAbsorb);

    pGEntity face;
    for (int i = 1; i <= 7; i++) {
        // Get the face
        face = GM_entityByTag(model, 2, i);

        // Add the face to the model association. Note that we passed
        // the Attribute Information Node into the Model Association
        // at the time when the Model Association was created. That prepares
        // the creation of the AttributeVoid on the face as soon as the
        // association process is started
        if (i == 6) {
            AMA_addGEntity(aSurf,face);
        } else if (i == 7) {
            AMA_addGEntity(aDynRup,face);
        } else {
            AMA_addGEntity(aAbsorb,face);
        }
    }

    // printFaceCoords(model);

    // ------------------------------ Set meshing parameters ------------------------------

    meshCase = MS_newMeshCase(model);

    // Set global mesh size
    pModelItem modelDomain = GM_domain(model);
    // ( <meshing case>, <entity>, <1=absolute, 2=relative>, <size>, <size expression> )
    MS_setMeshSize(meshCase, modelDomain, 1, 5000, NULL);

    // Set mesh size on fault
    face = GM_entityByTag(model, 2, 7);
    MS_setMeshSize(meshCase, face, 1, 400, NULL);

    // Set gradation relative
    MS_setGlobalSizeGradationRate(meshCase, 0.2);

    // Set target skewness
    MS_setVolumeShapeMetric(meshCase, modelDomain, ShapeMetricType_Skewness, 0.75);
}

void SimModSuite::loadSTL(const char *filename){
    pMesh mesh = M_new(0,0);
    pDiscreteModel d_model = 0;
    if(M_importFromSTLFile(mesh, filename, 0L)) { //check for error
        logError() << "Error importing file";
        M_release(mesh);
        return;
    }

    // check the input mesh for intersections
    // this call must occur before the discrete model is created
    if(MS_checkMeshIntersections(mesh, 0, 0L)) {
        logError() << "There are intersections in the input mesh";
        M_release(mesh);
        return;
    }

    // create the Discrete model
    d_model = DM_createFromMesh(mesh, 1, 0L);
    if(!d_model) { //check for error
        logError() << "Error creating Discrete model from mesh";
        M_release(mesh);
        return;
    }

    // define the Discrete model
    DM_findEdgesByFaceNormalsDegrees(d_model, 70, 0L);
    DM_eliminateDanglingEdges(d_model, 0L);
    if(DM_completeTopology(d_model, 0L)) { //check for error
        logError() << "Error completing Discrete model topology";
        M_release(mesh);
        GM_release(d_model);
        return;
    }

    // Print out information about the model
    logInfo(PMU_rank()) << "Number of model vertices: " << GM_numVertices(d_model);
    logInfo(PMU_rank()) << "Number of model edges: " << GM_numEdges(d_model);
    logInfo(PMU_rank()) << "Number of model faces: " << GM_numFaces(d_model);
    logInfo(PMU_rank()) << "Number of model regions: " << GM_numRegions(d_model);

    m_model = d_model;

    // Since we told the Discrete model to use the input mesh, we release our
    // pointer to it.  It will be fully released when the Discrete model is released.
    M_release(mesh);
}

void SimModSuite::loadCAD(const char* modFile, const char* cadFile){
    // Load CAD
    std::string sCadFile;
    if (cadFile)
        sCadFile = cadFile;
    else {
        sCadFile = modFile;
        utils::StringUtils::replaceLast(sCadFile, ".smd", "_nat.x_t");
    }
    pNativeModel nativeModel = 0L;
    if (utils::Path(sCadFile).exists())
        nativeModel = ParasolidNM_createFromFile(sCadFile.c_str(), 0);

    m_model = GM_load(modFile, nativeModel, 0L);

    if (nativeModel)
        NM_release(nativeModel);

    // check for model errors
    pPList modelErrors = PList_new();
    if (!GM_isValid(m_model, 0, modelErrors))
            // TODO print more detail about errors
            logError() << "Input model is not valid";
    PList_delete(modelErrors);
}
