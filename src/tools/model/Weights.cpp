#include "Weights.h"

#include <limits>
#include <maMesh.h>
#include <utils/logger.h>
#include <parmetis.h>
#include <fstream>
#include <cstring>
#include <cassert>

#include "SeismicVelocity.h"

unsigned getCluster(double timestep, double globalMinTimestep, unsigned rate)
{
  if (rate == 1) {
    return 0;
  }
  
  double upper;
  upper = rate * globalMinTimestep;
  
  unsigned cluster = 0;
  while (upper <= timestep) {
    upper *= rate;
    ++cluster;
  }
  return cluster;
}

int ipow(int x, int y) {
  assert(y > 0);

  if (y == 0) {
    return 1;
  }
  int result = x;
  while(--y) {
    result *= x;
  }
  return result;
}

idx_t* computeVertexWeights(apf::Mesh2* mesh, idx_t& ncon, int timestepRate, int drToCellRatio, bool enableDRWeights, char const* velocityModel)
{
  unsigned nLocalElements = apf::countOwned(mesh, 3);
  double* timesteps = new double[nLocalElements];
  int* dynamicRupture = new int[nLocalElements];
  double localMinTimestep = std::numeric_limits<double>::max();
  double localMaxTimestep = std::numeric_limits<double>::min();
  double drLocalMinTimestep = std::numeric_limits<double>::max();
  double globalMinTimestep, globalMaxTimestep, drGlobalMinTimestep;
  int localNumDrFaces = 0, globalNumDrFaces;
  
  std::fill(dynamicRupture, dynamicRupture + nLocalElements, 0);

  if (strlen(velocityModel) > 0) {
    double (*pWaveVelocityFun)(int,double,double,double);
    if (strcmp(velocityModel, "landers61") == 0) {
      pWaveVelocityFun = &landers61;
    } else if (strcmp(velocityModel, "sumatra1223") == 0) {
      pWaveVelocityFun = &sumatra1223;
    } else if (strcmp(velocityModel, "sumatra1224") == 0) {
      pWaveVelocityFun = &sumatra1224;
    } else {
      std::cerr << "Error: Unknown velocity model." << std::endl;
      exit(-1);
    }

    // Compute barycenter of each tetrahedron
    unsigned iElem = 0;
    apf::MeshTag* groupTag = mesh->findTag("group");
    apf::MeshIterator* it = mesh->begin(3);
    while (apf::MeshEntity* element = mesh->iterate(it)) {
      apf::Downward vertices;
      mesh->getDownward(element, 0, vertices);
      apf::Vector3 barycenter(0.,0.,0.);
      for (unsigned v = 0; v < 4; ++v) {
        apf::Vector3 x;
        mesh->getPoint(vertices[v], 0, x);
        barycenter += x * 0.25;
      }
      int group = -1;
      if (groupTag && mesh->hasTag(element, groupTag)) {
        mesh->getIntTag(element, groupTag, &group);
      }
      timesteps[iElem] = pWaveVelocityFun(group, barycenter.x(), barycenter.y(), barycenter.z());
      if (timesteps[iElem] < 0.0) {
        std::cerr << "Negative p wave velocity encountered." << std::endl;
        exit(-1);
      }
      ++iElem;
    }
    mesh->end(it);
  } else {
    std::fill(timesteps, timesteps + nLocalElements, 1.0);
  }
  
  // Compute timesteps
  unsigned iElem = 0;
  apf::MeshIterator* it = mesh->begin(3);
	while (apf::MeshEntity* element = mesh->iterate(it)) {
    timesteps[iElem] = ma::getInsphere(mesh, element) / timesteps[iElem];

    localMinTimestep = std::min(localMinTimestep, timesteps[iElem]);
    localMaxTimestep = std::max(localMaxTimestep, timesteps[iElem]);
    ++iElem;
	}
	mesh->end(it);
	MPI_Allreduce(&localMinTimestep, &globalMinTimestep, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(&localMaxTimestep, &globalMaxTimestep, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  
  apf::MeshTag* boundaryTag = mesh->findTag("boundary condition");
  it = mesh->begin(3);
  iElem = 0;
	while (apf::MeshEntity* element = mesh->iterate(it)) {
    apf::Downward faces;
    mesh->getDownward(element, 2, faces);
    for (unsigned f = 0; f < 4; ++f) {
      if (boundaryTag && mesh->hasTag(faces[f], boundaryTag)) {
        int boundary;
        mesh->getIntTag(faces[f], boundaryTag, &boundary);
        dynamicRupture[iElem] += (boundary == 3) ? 1 : 0;
      }
    }
    if (dynamicRupture[iElem] > 0) {
      drLocalMinTimestep = std::min(drLocalMinTimestep, timesteps[iElem]);
      localNumDrFaces = 1;
    } else {
      localNumDrFaces = 0;
    }
    ++iElem;
  }
	mesh->end(it);
  
  MPI_Allreduce(&drLocalMinTimestep, &drGlobalMinTimestep, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);  
  
  if (enableDRWeights) {
    MPI_Allreduce(&localNumDrFaces, &globalNumDrFaces, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);  
    ncon = (globalNumDrFaces > 0) ? 2 : 1;
  } else {
    ncon = 1;
  }
  unsigned maxCluster = getCluster(globalMaxTimestep, globalMinTimestep, timestepRate);
  idx_t* vwgt = new idx_t[ncon*nLocalElements];
  for (iElem = 0; iElem < nLocalElements; ++iElem) {
    double timestep = (dynamicRupture[iElem] == 0) ? timesteps[iElem] : drGlobalMinTimestep;
    // Actually the plus cell does all the work but I think this cannot
    // be adequately modeled here.
    vwgt[ncon*iElem] = (1 + drToCellRatio*dynamicRupture[iElem]) * ipow(timestepRate, maxCluster - getCluster(timestep, globalMinTimestep, timestepRate));
  }
  if (ncon > 1) {
    for (iElem = 0; iElem < nLocalElements; ++iElem) {
      vwgt[ncon*iElem+1] = (dynamicRupture[iElem] > 0) ? 1 : 0;
    }
  }
  
  delete[] timesteps;
  delete[] dynamicRupture;
  
  return vwgt;
}

idx_t* computeEdgeWeights(apf::Mesh2* mesh, int const* dualGraph, idx_t nEdges)
{
  unsigned nLocalElements = apf::countOwned(mesh, 3);
  idx_t* adjwgt = new idx_t[nEdges];
  std::fill(adjwgt, adjwgt + nEdges, 1.0);
  apf::MeshTag* boundaryTag = mesh->findTag("boundary condition");
	unsigned int pos = 0;
	apf::MeshIterator* it = mesh->begin(3);
	for (unsigned int e = 0; e < nLocalElements; ++e) {
    apf::MeshEntity* element = mesh->iterate(it);
    apf::Downward faces;
    mesh->getDownward(element, 2, faces);
    for (unsigned int f = 0; f < 4; ++f) {
      if (dualGraph[e*4 + f] >= 0) {
        if (boundaryTag && mesh->hasTag(faces[f], boundaryTag)) {
          int boundary;
          mesh->getIntTag(faces[f], boundaryTag, &boundary);
          if (boundary == 3) {
            adjwgt[pos] = 100.0;
          }
        }
        ++pos;
      }
    }
	}
  mesh->end(it);
  return adjwgt;
}
