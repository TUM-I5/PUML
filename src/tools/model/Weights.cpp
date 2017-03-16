#include "Weights.h"

#include <limits>
#include <maMesh.h>
#include <PCU.h>
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

void computeTimesteps(apf::Mesh2* mesh, char const* velocityModel, double& globalMinTimestep, double& globalMaxTimestep)
{
  double localMinTimestep = std::numeric_limits<double>::max();
  double localMaxTimestep = std::numeric_limits<double>::min();
	apf::MeshTag* timestepTag = mesh->createDoubleTag("timestep", 1);

  if (strlen(velocityModel) > 0) {
    double (*pWaveVelocityFun)(int,double,double,double);
    if (strcmp(velocityModel, "landers61") == 0) {
      pWaveVelocityFun = &landers61;
    } else if (strcmp(velocityModel, "sumatra1223") == 0) {
	pWaveVelocityFun = 0L;
	logError() << "Obsolete velocity model, use \"sumatra1223_high\" or \"sumatra1223_low\"";
    } else if (strcmp(velocityModel, "sumatra1223_high") == 0) {
      pWaveVelocityFun = &sumatra1223_high;
    } else if (strcmp(velocityModel, "sumatra1223_low") == 0) {
      pWaveVelocityFun = &sumatra1223_low;
    } else if (strcmp(velocityModel, "sumatra1224") == 0) {
      pWaveVelocityFun = &sumatra1224;
    } else {
      std::cerr << "Error: Unknown velocity model." << std::endl;
      exit(-1);
    }

    // Compute barycenter of each tetrahedron
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
      double timestep = pWaveVelocityFun(group, barycenter.x(), barycenter.y(), barycenter.z());
      if (timestep < 0.0) {
        std::cerr << "Negative p wave velocity encountered." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
      }
      mesh->setDoubleTag(element, timestepTag, &timestep);
    }
    mesh->end(it);
  } else {
    apf::MeshIterator* it = mesh->begin(3);
    double one = 1.0;
    while (apf::MeshEntity* element = mesh->iterate(it)) {
      mesh->setDoubleTag(element, timestepTag, &one);
    }
    mesh->end(it);
  }

  apf::MeshIterator* it = mesh->begin(3);
	while (apf::MeshEntity* element = mesh->iterate(it)) {
    double timestep;
    mesh->getDoubleTag(element, timestepTag, &timestep);
    timestep = ma::getInsphere(mesh, element) / timestep;
    mesh->setDoubleTag(element, timestepTag, &timestep);

    localMinTimestep = std::min(localMinTimestep, timestep);
    localMaxTimestep = std::max(localMaxTimestep, timestep);
	}
	mesh->end(it);
	MPI_Allreduce(&localMinTimestep, &globalMinTimestep, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(&localMaxTimestep, &globalMaxTimestep, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
}

void countDynamicRuptureFaces(apf::Mesh2* mesh, int& globalNumDrFaces)
{
  int localNumDrFaces = 0;
	apf::MeshTag* dynRupTag = mesh->createIntTag("dynamicRupture", 1);
  apf::MeshTag* boundaryTag = mesh->findTag("boundary condition");

  apf::MeshIterator* it = mesh->begin(3);
	while (apf::MeshEntity* element = mesh->iterate(it)) {
    apf::Downward faces;
    mesh->getDownward(element, 2, faces);
    int dynamicRupture = 0;
    for (unsigned f = 0; f < 4; ++f) {
      if (boundaryTag && mesh->hasTag(faces[f], boundaryTag)) {
        int boundary;
        mesh->getIntTag(faces[f], boundaryTag, &boundary);
        dynamicRupture += (boundary == 3) ? 1 : 0;
      }
    }
    if (dynamicRupture > 0) {
      ++localNumDrFaces;
    }
    mesh->setIntTag(element, dynRupTag, &dynamicRupture);
  }
	mesh->end(it);

  MPI_Allreduce(&localNumDrFaces, &globalNumDrFaces, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
}

int enforceDynamicRuptureGTS(apf::Mesh2* mesh)
{
  int numberOfReductions = 0;

	apf::MeshTag* dynRupTag = mesh->findTag("dynamicRupture");
	apf::MeshTag* clusterTag = mesh->findTag("timeCluster");

  int localMinCluster = std::numeric_limits<int>::max(), globalMinCluster;
  apf::MeshIterator* it = mesh->begin(3);
  while (apf::MeshEntity* element = mesh->iterate(it)) {
    int dynamicRupture;
    mesh->getIntTag(element, dynRupTag, &dynamicRupture);
    if (dynamicRupture > 0) {
      int timeCluster;
      mesh->getIntTag(element, clusterTag, &timeCluster);
      localMinCluster = std::min(localMinCluster, timeCluster);
    }
  }
	mesh->end(it);

  MPI_Allreduce(&localMinCluster, &globalMinCluster, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  it = mesh->begin(3);
  while (apf::MeshEntity* element = mesh->iterate(it)) {
    int dynamicRupture, timeCluster;
    mesh->getIntTag(element, dynRupTag, &dynamicRupture);
    mesh->getIntTag(element, clusterTag, &timeCluster);
    if (dynamicRupture > 0 && timeCluster != globalMinCluster) {
      mesh->setIntTag(element, clusterTag, &globalMinCluster);
      ++numberOfReductions;
    }
  }
	mesh->end(it);

  return numberOfReductions;
}

int normalizeClustering(apf::Mesh2* mesh, int maxDifference = 1)
{
	apf::MeshTag* clusterTag = mesh->findTag("timeCluster");
  apf::MeshTag* boundaryTag = mesh->findTag("boundary condition");

  int numberOfReductions = 0;

  apf::MeshIterator* it = mesh->begin(3);
  while (apf::MeshEntity* element = mesh->iterate(it)) {
    int timeCluster;
    mesh->getIntTag(element, clusterTag, &timeCluster);

    apf::Downward faces;
    mesh->getDownward(element, 2, faces);
    for (unsigned f = 0; f < 4; ++f) {
      int boundary = -1;
      if (mesh->hasTag(faces[f], boundaryTag)) {
        mesh->getIntTag(faces[f], boundaryTag, &boundary);
      }
      // Continue for regular, dynamic rupture, and periodic boundary cells
      if (boundary == -1 || boundary == 3 || boundary == 6) {
        // We treat MPI neighbours later
        if (!mesh->isShared(faces[f])) {
          apf::Up elements;
          mesh->getUp(faces[f], elements);

          if (elements.n != 2) {
            std::cerr << "Could not find a face neighbour." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
          }

          apf::MeshEntity* neighbour = (elements.e[0] == element) ? elements.e[1] : elements.e[0];
          int otherTimeCluster;
          mesh->getIntTag(neighbour, clusterTag, &otherTimeCluster);

          if (timeCluster > otherTimeCluster + maxDifference) {
            timeCluster = otherTimeCluster + maxDifference;
            ++numberOfReductions;
          }
        }
      }
    }
    mesh->setIntTag(element, clusterTag, &timeCluster);
  }
  mesh->end(it);

  PCU_Comm_Begin();
	it = mesh->begin(3);
	while (apf::MeshEntity* element = mesh->iterate(it)) {
		apf::Downward faces;
		mesh->getDownward(element, 2, faces);

		for (unsigned int i = 0; i < 4; i++) {
			if (!mesh->isShared(faces[i])) {
				continue;
      }
			apf::Copy other = apf::getOtherCopy(mesh, faces[i]);
			PCU_COMM_PACK(other.peer, other.entity);
			int timeCluster;
			mesh->getIntTag(element, clusterTag, &timeCluster);
			PCU_COMM_PACK(other.peer, timeCluster);
		}
	}
	mesh->end(it);
	PCU_Comm_Send();

	while (PCU_Comm_Receive()) {
		apf::MeshEntity* face;
		PCU_COMM_UNPACK(face);
		int otherTimeCluster, timeCluster;
		PCU_Comm_Unpack(&otherTimeCluster, sizeof(otherTimeCluster));

    apf::Up elements;
    mesh->getUp(face, elements);
    if (elements.n != 1) {
      std::cerr << "That is unexpected." << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
    mesh->getIntTag(elements.e[0], clusterTag, &timeCluster);
    if (timeCluster > otherTimeCluster + maxDifference) {
      timeCluster = otherTimeCluster + maxDifference;
      mesh->setIntTag(elements.e[0], clusterTag, &timeCluster);
      ++numberOfReductions;
    }
	}

  return numberOfReductions;
}

idx_t* computeVertexWeights(apf::Mesh2* mesh, idx_t& ncon, int timestepRate, int drToCellRatio, bool enableDRWeights, char const* velocityModel)
{
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  unsigned nLocalElements = apf::countOwned(mesh, 3);
  idx_t* vwgt = new idx_t[ncon*nLocalElements];
  double globalMinTimestep, globalMaxTimestep;
  int globalNumDrFaces;

  computeTimesteps(mesh, velocityModel, globalMinTimestep, globalMaxTimestep);
  countDynamicRuptureFaces(mesh, globalNumDrFaces);

	apf::MeshTag* dynRupTag = mesh->findTag("dynamicRupture");
	apf::MeshTag* timestepTag = mesh->findTag("timestep");
	apf::MeshTag* clusterTag = mesh->createIntTag("timeCluster", 1);
  apf::MeshIterator* it = mesh->begin(3);
  while (apf::MeshEntity* element = mesh->iterate(it)) {
    double timestep;
    mesh->getDoubleTag(element, timestepTag, &timestep);
    int timeCluster = getCluster(timestep, globalMinTimestep, timestepRate);
    mesh->setIntTag(element, clusterTag, &timeCluster);
  }
  mesh->end(it);

  if (timestepRate > 1) {
    int totalNumberOfReductions = 0;
    int localNumberOfReductions, globalNumberOfReductions;
    do {
      localNumberOfReductions = 0;
      if (globalNumDrFaces > 0) {
        localNumberOfReductions += enforceDynamicRuptureGTS(mesh);
      }
      localNumberOfReductions += normalizeClustering(mesh);

      MPI_Allreduce(&localNumberOfReductions, &globalNumberOfReductions, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      totalNumberOfReductions += globalNumberOfReductions;
    } while (globalNumberOfReductions > 0);

    logInfo(rank) << "Number of reductions:" << totalNumberOfReductions;
  }

  if (enableDRWeights) {
    ncon = (globalNumDrFaces > 0) ? 2 : 1;
  } else {
    ncon = 1;
  }

  unsigned maxCluster = getCluster(globalMaxTimestep, globalMinTimestep, timestepRate);
  unsigned iElem = 0;
  it = mesh->begin(3);
  while (apf::MeshEntity* element = mesh->iterate(it)) {
    int timeCluster, dynamicRupture;
    mesh->getIntTag(element, clusterTag, &timeCluster);
    mesh->getIntTag(element, dynRupTag, &dynamicRupture);
    if (ncon > 1) {
      vwgt[ncon*iElem] = ipow(timestepRate, maxCluster - timeCluster);
      vwgt[ncon*iElem+1] = (dynamicRupture > 0) ? 1 : 0;
    } else {
      // Actually the plus cell does all the work but I think this cannot
      // be adequately modeled here.
      vwgt[ncon*iElem] = (1 + drToCellRatio*dynamicRupture) * ipow(timestepRate, maxCluster - timeCluster);
    }
    ++iElem;
  }
	mesh->end(it);

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
