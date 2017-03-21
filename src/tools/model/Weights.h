#ifndef WEIGHTS_H_
#define WEIGHTS_H_

#include <apfMesh2.h>
#include <parmetis.h>

unsigned getCluster(double timestep, double globalMinTimestep, unsigned rate);
idx_t* computeVertexWeights(apf::Mesh2* mesh, idx_t& ncon, int timestepRate, int drToCellRatio, bool enableDRWeights, char const* velocityModel, unsigned& maxCluster);
idx_t* computeEdgeWeights(apf::Mesh2* mesh, int const* dualGraph, idx_t nEdges, int timestepRate, unsigned maxCluster);

#endif
