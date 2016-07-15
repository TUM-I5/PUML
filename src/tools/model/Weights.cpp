#include "Weights.h"

#include <limits>
#include <maMesh.h>
#include <utils/logger.h>
#include <projects.h>
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

idx_t* computeVertexWeights(apf::Mesh2* mesh, char const* sourceCoordSystem)
{
  unsigned nLocalElements = apf::countOwned(mesh, 3);
  double* lat = new double[nLocalElements];
  double* lon = new double[nLocalElements];
  double* height = new double[nLocalElements];
  double* timesteps = new double[nLocalElements];
  bool* dynamicRupture = new bool[nLocalElements];
  double localMinTimestep = std::numeric_limits<double>::max();
  double localMaxTimestep = std::numeric_limits<double>::min();
  double globalMinTimestep, globalMaxTimestep;
  
  std::fill(dynamicRupture, dynamicRupture + nLocalElements, false);

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
  }
	mesh->end(it);
  
  projPJ pj_lonlat;
  projPJ pj_mesh;
  // Transform from mesh coordinate system to latitude, longitude, height
  if (!(pj_mesh = pj_init_plus(sourceCoordSystem))) { 
    projPrintLastError();
  }
  if (!(pj_lonlat = pj_init_plus("+proj=lonlat +datum=WGS84 +units=m +no_defs"))) { 
    projPrintLastError();
  }
  pj_transform(pj_mesh, pj_lonlat, nLocalElements, 1, lat, lon, height);

  pj_free(pj_lonlat);
  pj_free(pj_mesh);
  
  // Compute maximum wave velocity (= P wave)
  get_material_parameters(lat, lon, height, nLocalElements, Vp, timesteps);
  
  delete[] lat;
  delete[] lon;
  delete[] height;
  
  // Compute timesteps
  iElem = 0;
  it = mesh->begin(3);
	while (apf::MeshEntity* element = mesh->iterate(it)) {
    timesteps[iElem] = ma::getInsphere(mesh, element) / timesteps[iElem]; // TODO: divide by wavespeed

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
        dynamicRupture[iElem] = dynamicRupture[iElem] || (boundary == 3);
      }
    }
    ++iElem;
  }
	mesh->end(it);
  
  unsigned maxCluster = getCluster(globalMaxTimestep, globalMinTimestep, 2);
  idx_t* vwgt = new idx_t[2*nLocalElements];
  for (iElem = 0; iElem < nLocalElements; ++iElem) {
    vwgt[2*iElem] = (1 << (maxCluster - getCluster(timesteps[iElem], globalMinTimestep, 2))); // Valid for rate 2
    vwgt[2*iElem+1] = (dynamicRupture[iElem]) ? 1 : 0;
  }
  
  delete[] timesteps;
  delete[] dynamicRupture;
  
  return vwgt;
}
