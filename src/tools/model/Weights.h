#ifndef WEIGHTS_H_
#define WEIGHTS_H_

#include <apfMesh2.h>
#include <parmetis.h>

unsigned getCluster(double timestep, double globalMinTimestep, unsigned rate);
idx_t* computeVertexWeights(apf::Mesh2* mesh, char const* sourceCoordSystem);

#endif
