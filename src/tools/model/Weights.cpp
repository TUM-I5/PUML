#include "Weights.h"

#include <limits>
#include <maMesh.h>
#include <utils/logger.h>
#include <projects.h>
#include <parmetis.h>
#include "SeismicVelocity.h"

void projPrintLastError()
{
  if (pj_errno != 0) {
    logError() << "Proj.4 error:" << pj_strerrno(pj_errno);
  }
}

unsigned getCluster(double timestep, double globalMinTimestep, unsigned rate)
{
  double upper;
  upper = rate * globalMinTimestep;
  
  unsigned cluster = 0;
  while (upper <= timestep) {
    upper *= rate;
    ++cluster;
  }
  return cluster;
}

idx_t* computeVertexWeights(apf::Mesh2* mesh, char const* sourceCoordSystem, idx_t& ncon, int drToCellRatio, bool enableDRWeights)
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
  
  if (strlen(sourceCoordSystem) > 0) {
    double* lat = new double[nLocalElements];
    double* lon = new double[nLocalElements];
    double* height = new double[nLocalElements];
  
    // Compute barycenter of each tetrahedron
    unsigned iElem = 0;
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
      lat[iElem] = barycenter.x();
      lon[iElem] = barycenter.y();
      height[iElem] = barycenter.z();
      ++iElem;
    }
    mesh->end(it);
    
    projPJ pj_lonlat;
    projPJ pj_mesh;
    // Transform from mesh coordinate system to latitude, longitude, height
    if (!(pj_mesh = pj_init_plus(sourceCoordSystem))) { 
      projPrintLastError();
    }
    if (!(pj_lonlat = pj_init_plus("+proj=latlon +datum=WGS84 +units=m +no_defs"))) {
      projPrintLastError();
    }

    pj_transform(pj_mesh, pj_lonlat, nLocalElements, 1, lat, lon, height);
    projPrintLastError();

    pj_free(pj_lonlat);
    pj_free(pj_mesh);
    
    // Compute maximum wave velocity (= P wave)
    get_material_parameters(lat, lon, height, nLocalElements, Vp, timesteps);
    
    delete[] lat;
    delete[] lon;
    delete[] height;
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
  unsigned maxCluster = getCluster(globalMaxTimestep, globalMinTimestep, 2);
  idx_t* vwgt = new idx_t[ncon*nLocalElements];
  for (iElem = 0; iElem < nLocalElements; ++iElem) {
    double timestep = (dynamicRupture[iElem] == 0) ? timesteps[iElem] : drGlobalMinTimestep;
    // Actually the plus cell does all the work but I think this cannot
    // be adequately modeled here.
    vwgt[ncon*iElem] = (1 + drToCellRatio*dynamicRupture[iElem]) * (1 << (maxCluster - getCluster(timestep, globalMinTimestep, 2))); // Valid for rate 2
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
