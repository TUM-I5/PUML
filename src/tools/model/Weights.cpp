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

#include "netcdf.h"

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
  assert(y >= 0);

  if (y == 0) {
    return 1;
  }
  int result = x;
  while(--y) {
    result *= x;
  }
  return result;
}

void check_err(const int stat, const int line, const char *file) {
  if (stat != NC_NOERR) {
    logError() << "line" << line << "of" << file << ":" << nc_strerror(stat) << std::endl;
  }
}

void readAsagiDerivedNetcdfData(int &nx, int &ny, int &nz, double * &xNc, double * &yNc,double * &zNc, double * &VpNc)
{
   int stat;
   int ncFile, varid,dimid;
   size_t na,nb,nc;

   stat = (nc_open("AsagiDerivedVp.nc", NC_NETCDF4, &ncFile));
   check_err(stat,__LINE__,__FILE__);

   stat = (nc_inq_dimid(ncFile, "x", &dimid));
   check_err(stat,__LINE__,__FILE__);
   stat = (nc_inq_dimlen(ncFile, dimid, &na));
   check_err(stat,__LINE__,__FILE__);
   stat = (nc_inq_dimid(ncFile, "y", &dimid));
   check_err(stat,__LINE__,__FILE__);
   stat = (nc_inq_dimlen(ncFile, dimid, &nb));
   check_err(stat,__LINE__,__FILE__);
   stat = (nc_inq_dimid(ncFile, "z", &dimid));
   check_err(stat,__LINE__,__FILE__);
   stat = (nc_inq_dimlen(ncFile, dimid, &nc));
   check_err(stat,__LINE__,__FILE__);

   nx=static_cast<int>(na);
   ny=static_cast<int>(nb);
   nz=static_cast<int>(nc);


   xNc = new double [nx];
   yNc = new double [ny];
   zNc = new double [nz];
   VpNc = new double [nx*ny*nz];

   stat = (nc_inq_varid(ncFile, "x", &varid));
   check_err(stat,__LINE__,__FILE__);
   stat = (nc_get_var(ncFile, varid, &xNc[0]));
   check_err(stat,__LINE__,__FILE__);
   stat = (nc_inq_varid(ncFile, "y", &varid));
   check_err(stat,__LINE__,__FILE__);
   stat = (nc_get_var(ncFile, varid, &yNc[0]));
   check_err(stat,__LINE__,__FILE__);
   stat = (nc_inq_varid(ncFile, "z", &varid));
   check_err(stat,__LINE__,__FILE__);
   stat = (nc_get_var(ncFile, varid, &zNc[0]));
   check_err(stat,__LINE__,__FILE__);
   stat = (nc_inq_varid(ncFile, "Vp", &varid));
   check_err(stat,__LINE__,__FILE__);
   stat = (nc_get_var(ncFile, varid, &VpNc[0]));
   check_err(stat,__LINE__,__FILE__);
   
}


void computeTimesteps(apf::Mesh2* mesh, char const* velocityModel, double& globalMinTimestep, double& globalMaxTimestep)
{
  double localMinTimestep = std::numeric_limits<double>::max();
  double localMaxTimestep = std::numeric_limits<double>::min();
	apf::MeshTag* timestepTag = mesh->createDoubleTag("timestep", 1);
  int nx,ny,nz;
  double *xNc = NULL;
  double *yNc = NULL;
  double *zNc = NULL;
  double *VpNc = NULL;

  if (strlen(velocityModel) > 0) {
    double (*pWaveVelocityFun)(int,double,double,double);
    if (strcmp(velocityModel, "AsagiDerived") == 0) {
      readAsagiDerivedNetcdfData(nx,ny,nz,xNc,yNc,zNc,VpNc);
      //pWaveVelocityFun has a different list of arguments, so we cant use the function pointer here
    } else if (strcmp(velocityModel, "landers61") == 0) {
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
      double timestep;
      if (strcmp(velocityModel, "AsagiDerived") == 0) {
         timestep = AsagiDerived(group, barycenter.x(), barycenter.y(), barycenter.z(), nx,ny,nz,xNc,yNc,zNc,VpNc);
      } else {
         timestep = pWaveVelocityFun(group, barycenter.x(), barycenter.y(), barycenter.z());
      }
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

int enforceMaximumDifference(apf::Mesh2* mesh, int maxDifference = 1)
{
	apf::MeshTag* clusterTag = mesh->findTag("timeCluster");
  apf::MeshTag* boundaryTag = mesh->findTag("boundary condition");
  apf::MeshTag* dynRupTag = mesh->findTag("dynamicRupture");

  int numberOfReductions = 0;

  apf::MeshIterator* it = mesh->begin(3);
  while (apf::MeshEntity* element = mesh->iterate(it)) {
    int timeCluster;
    mesh->getIntTag(element, clusterTag, &timeCluster);

    apf::Downward faces;
    mesh->getDownward(element, 2, faces);
    for (unsigned f = 0; f < 4; ++f) {
      int difference = maxDifference;
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
            std::cerr << "For that element: boundary = " << boundary << std::endl;
            apf::Vector3 tmp;
            tmp = getLinearCentroid(mesh,elements.e[0]);
            std::cerr << "linear centroid"<< tmp[0]<< " " << tmp[1]<< " "<<tmp[2]<< std::endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
          }

          apf::MeshEntity* neighbour = (elements.e[0] == element) ? elements.e[1] : elements.e[0];
          int otherTimeCluster;
          mesh->getIntTag(neighbour, clusterTag, &otherTimeCluster);
          
          if (boundary == 3) {
            difference = 0;
          }

          if (timeCluster > otherTimeCluster + difference) {
            timeCluster = otherTimeCluster + difference;
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
    int difference = maxDifference;

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
    int boundary = -1;
    if (mesh->hasTag(face, boundaryTag)) {
      mesh->getIntTag(face, boundaryTag, &boundary);
    }
    if (boundary == -1 || boundary == 3 || boundary == 6) {
      mesh->getIntTag(elements.e[0], clusterTag, &timeCluster);
      if (boundary == 3) {
        difference = 0;
      }
      if (timeCluster > otherTimeCluster + difference) {
        timeCluster = otherTimeCluster + difference;
        mesh->setIntTag(elements.e[0], clusterTag, &timeCluster);
        ++numberOfReductions;
      }
    }
  }

  return numberOfReductions;
}

idx_t* computeVertexWeights(apf::Mesh2* mesh, idx_t& ncon, int timestepRate, int drToCellRatio, bool enableDRWeights, char const* velocityModel, unsigned& maxCluster)
{
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  unsigned nLocalElements = apf::countOwned(mesh, 3);
  double globalMinTimestep, globalMaxTimestep;
  int globalNumDrFaces;

  logInfo(rank) << "Computing timesteps";
  computeTimesteps(mesh, velocityModel, globalMinTimestep, globalMaxTimestep);
  countDynamicRuptureFaces(mesh, globalNumDrFaces);

  logInfo(rank) << "Determining time clusters";
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
    logInfo(rank) << "Normalizing time clusters";
    int totalNumberOfReductions = 0;
    int localNumberOfReductions, globalNumberOfReductions;
    do {
      localNumberOfReductions = enforceMaximumDifference(mesh);

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
  
  idx_t* vwgt = new idx_t[ncon*nLocalElements];
  maxCluster = getCluster(globalMaxTimestep, globalMinTimestep, timestepRate);
  int* localClusterHistogram = new int[maxCluster+1];
  int* globalClusterHistogram;
  if (rank == 0) {
    globalClusterHistogram = new int[maxCluster+1];
  }
  for (int cluster = 0; cluster <= maxCluster; ++cluster) {
    localClusterHistogram[cluster] = 0;
  }
  unsigned iElem = 0;
  it = mesh->begin(3);
  while (apf::MeshEntity* element = mesh->iterate(it)) {
    int timeCluster, dynamicRupture;
    mesh->getIntTag(element, clusterTag, &timeCluster);
    mesh->getIntTag(element, dynRupTag, &dynamicRupture);
    ++localClusterHistogram[timeCluster];
    
    // Actually the plus cell does all the work but I think this cannot
    // be adequately modeled here.
    vwgt[ncon*iElem] = (1 + drToCellRatio*dynamicRupture) * ipow(timestepRate, maxCluster - timeCluster);
    if (ncon > 1) {
      vwgt[ncon*iElem+1] = (dynamicRupture > 0) ? 1 : 0;
    }
    ++iElem;
  }
	mesh->end(it);
  
  MPI_Reduce(localClusterHistogram, globalClusterHistogram, maxCluster+1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (rank == 0) {
    logInfo(rank) << "Number of elements in time clusters:";
    for (int cluster = 0; cluster <= maxCluster; ++cluster) {
      logInfo(rank) << utils::nospace << cluster << ":" << utils::space << globalClusterHistogram[cluster];
    }
    delete[] globalClusterHistogram;
  }  
  delete[] localClusterHistogram;

  return vwgt;
}

void determineEdgeWeights(apf::Mesh2* mesh, int timestepRate, unsigned maxCluster)
{
	apf::MeshTag* clusterTag = mesh->findTag("timeCluster");
  apf::MeshTag* boundaryTag = mesh->findTag("boundary condition");
  apf::MeshTag* edgeWeightTag = mesh->findTag("edgeWeight");

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
            std::cerr << "For that element: boundary = " << boundary << std::endl;
            apf::Vector3 tmp;
            tmp = getLinearCentroid(mesh,elements.e[0]);
            std::cerr << "linear centroid"<< tmp[0]<< " " << tmp[1]<< " "<<tmp[2]<< std::endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
          }

          apf::MeshEntity* neighbour = (elements.e[0] == element) ? elements.e[1] : elements.e[0];
          int otherTimeCluster;
          mesh->getIntTag(neighbour, clusterTag, &otherTimeCluster);

          int edgeWeight = ipow(timestepRate, maxCluster - std::max(timeCluster, otherTimeCluster));
          mesh->setIntTag(faces[f], edgeWeightTag, &edgeWeight);
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
    int boundary = -1;
    if (mesh->hasTag(face, boundaryTag)) {
      mesh->getIntTag(face, boundaryTag, &boundary);
    }
    if (boundary == -1 || boundary == 3 || boundary == 6) {
      mesh->getIntTag(elements.e[0], clusterTag, &timeCluster);
      int edgeWeight = ipow(timestepRate, maxCluster - std::max(timeCluster, otherTimeCluster));
      mesh->setIntTag(face, edgeWeightTag, &edgeWeight);
    }
  }
}

idx_t* computeEdgeWeights(apf::Mesh2* mesh, int const* dualGraph, idx_t nEdges, int timestepRate, unsigned maxCluster)
{
	apf::MeshTag* edgeWeightTag = mesh->createIntTag("edgeWeight", 1);
  determineEdgeWeights(mesh, timestepRate, maxCluster);
  
  unsigned nLocalElements = apf::countOwned(mesh, 3);
  idx_t* adjwgt = new idx_t[nEdges];
  std::fill(adjwgt, adjwgt + nEdges, 1.0);
	unsigned int pos = 0;
	apf::MeshIterator* it = mesh->begin(3);
	for (unsigned int e = 0; e < nLocalElements; ++e) {
    apf::MeshEntity* element = mesh->iterate(it);
    apf::Downward faces;
    mesh->getDownward(element, 2, faces);
    for (unsigned int f = 0; f < 4; ++f) {
      if (dualGraph[e*4 + f] >= 0) {
        if (!mesh->hasTag(faces[f], edgeWeightTag)) {
          std::cerr << "Something went terribly wrong." << std::endl;
          MPI_Abort(MPI_COMM_WORLD, -1);
        }
        int edgeWeight;
        mesh->getIntTag(faces[f], edgeWeightTag, &edgeWeight);
        adjwgt[pos] = edgeWeight;
        ++pos;
      }
    }
	}
  mesh->end(it);
  return adjwgt;
}
