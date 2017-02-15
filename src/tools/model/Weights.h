#ifndef WEIGHTS_H_
#define WEIGHTS_H_

#include <apfMesh2.h>
#include <parmetis.h>

unsigned getCluster(double timestep, double globalMinTimestep, unsigned rate);
idx_t* computeVertexWeights(apf::Mesh2* mesh, char const* sourceCoordSystem, idx_t& ncon, int drToCellRatio, bool enableDRWeights);
idx_t* computeEdgeWeights(apf::Mesh2* mesh, int const* dualGraph, idx_t nEdges);

#endif
