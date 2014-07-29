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

#ifndef H_SimModelerUtil
#define H_SimModelerUtil

#ifndef AdvMeshing_EXPORT
#define AdvMeshing_EXPORT
#endif

#include "SimAttribute.h"
#include "ModelTypes.h"

/* The fields in this sturcture are set mainly based on the "Surface Meshing" and "Volume Meshing"
   attributes in the meshing case. */
struct AdvMeshing_EXPORT MeshingOptions {
  bool surfaceRun;                // whether to run surface meshing
  bool surfaceDoFixIntersections; // whether to run fix self intersections
  bool surfaceDoCurve;            // whether to curve surface mesh

  // These correspond as indicated to the fields in the "Surface Meshing" attribute
  int surfaceSmoothingLevel;      // "Smoothing Level"
  int surfaceSmoothingType;       // "Smoothing Type"
  int surfaceFixIntersections;    // "Fix Intersections"
  int surfaceOptimization;        // "Optimization"
  int surfaceEnforceGradation;    // "Enforce Spatial Gradation"
  double surfaceFaceRotationLimit;// "Discrete Face Rotation Limit"
  int surfaceCurveType;           // ignore

  bool volumeRun;                 // whether to run volume meshing
  bool volumeDoStructured;        // always true
  bool volumeDoUnstructured;      // always true
  bool volumeDoCurve;             // whether to run mesh curving

  // These correspond as indicated to the fields in the "Volume Meshing" attribute
  int volumeEnforceSize;          // "Enforce Size"
  int volumeSmoothingLevel;       // "Smoothing Level"
  int volumeSmoothingType;        // "Smoothing Type"
  int volumeOptimization;         // "Optimization"
  int volumeCurveType;            // ignore
};

AdvMeshing_EXPORT void MS_setupSimModelerMeshCase(pACase sourceCase, pACase meshCase, MeshingOptions *options);

#endif
