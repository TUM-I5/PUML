/**
 * @file
 *  This file is part of PUML
 *
 *  For conditions of distribution and use, please see the copyright
 *  notice in the file 'COPYING' at the root directory of this package
 *  and the copyright notice at https://github.com/TUM-I5/PUML
 *
 * @copyright 2015 Technische Universitaet Muenchen
 * @author Sebastian Rettenberger <rettenbs@in.tum.de>
 */

#ifndef NETCDF_MESH_H
#define NETCDF_MESH_H

#ifdef PARALLEL
#include <mpi.h>
#endif // PARALLEL

#include <netcdf.h>
#ifdef PARALLEL
#include <netcdf_par.h>
#endif // PARALLEL

#include <apfConvert.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <gmi_null.h>

#include "utils/logger.h"

#include "MeshInput.h"
#include "NetCDFPartition.h"

/**
 * Read PUMGen generated mesh files
 */
class NetCDFMesh : public MeshInput
{
public:
	NetCDFMesh(const char* meshFile, MPI_Comm comm = MPI_COMM_WORLD)
	{
		int rank = 0;
		int nProcs = 1;
#ifdef PARALLEL
		MPI_Comm_rank(comm, &rank);
		MPI_Comm_size(comm, &nProcs);
#endif // PARALLEL

		gmi_register_null();
		gmi_model* model = gmi_load(".null");
		m_mesh = apf::makeEmptyMdsMesh(model, 3, false);

		int ncFile;
#ifdef PARALLEL
		checkNcError(nc_open_par(meshFile, NC_MPIIO,
				comm, MPI_INFO_NULL, &ncFile));
#else // PARALLEL
		checkNcError(nc_open(meshFile, 0, &ncFile));
#endif // PARALLEL

		// Get number of partitions
		int ncDimPart;
		checkNcError(nc_inq_dimid(ncFile, "partitions", &ncDimPart));
		size_t nPartitions;
		checkNcError(nc_inq_dimlen(ncFile, ncDimPart, &nPartitions));

		// Local partitions
		unsigned int nMaxLocalPart = (nPartitions + nProcs - 1) / nProcs;
		unsigned int nLocalPart = nMaxLocalPart;
		if (nPartitions < (rank+1) * nMaxLocalPart)
			nLocalPart = std::max(0, static_cast<int>(nPartitions - rank * nMaxLocalPart));

#ifdef PARALLEL
		MPI_Comm commIO;
		MPI_Comm_split(MPI_COMM_WORLD, (nLocalPart > 0 ? 0 : MPI_UNDEFINED), 0, &commIO);

		// Reopen netCDF file with correct communicator
		checkNcError(nc_close(ncFile));

		if (nLocalPart > 0)
			checkNcError(nc_open_par(meshFile, NC_MPIIO,
					commIO, MPI_INFO_NULL, &ncFile));

		PCU_Switch_Comm(commIO);
#endif // PARALLEL

		unsigned int nLocalElements = 0;
		unsigned int nLocalVertices = 0;
		int* elements = 0L;
		double* vertices = 0L;
		int* boundaries = 0L;
		int* groups = 0L;

		if (nLocalPart > 0) {
			// Create netCDF variables
			int ncVarElemSize;
			checkNcError(nc_inq_varid(ncFile, "element_size", &ncVarElemSize));
			collectiveAccess(ncFile, ncVarElemSize);

			int ncVarElemVertices;
			checkNcError(nc_inq_varid(ncFile, "element_vertices", &ncVarElemVertices));
			collectiveAccess(ncFile, ncVarElemVertices);

			int ncVarElemNeighborRanks;
			checkNcError(nc_inq_varid(ncFile, "element_neighbor_ranks", &ncVarElemNeighborRanks));
			collectiveAccess(ncFile, ncVarElemNeighborRanks);

			int ncVarElemSideOrientations;
			checkNcError(nc_inq_varid(ncFile, "element_side_orientations", &ncVarElemSideOrientations));
			collectiveAccess(ncFile, ncVarElemSideOrientations);

			int ncVarElemMPIIndices;
			checkNcError(nc_inq_varid(ncFile, "element_mpi_indices", &ncVarElemMPIIndices));
			collectiveAccess(ncFile, ncVarElemMPIIndices);

			int ncVarElemBoundaries;
			checkNcError(nc_inq_varid(ncFile, "element_boundaries", &ncVarElemBoundaries));
			collectiveAccess(ncFile, ncVarElemBoundaries);

			int ncVarElemGroup;
			checkNcError(nc_inq_varid(ncFile, "element_group", &ncVarElemGroup));
			collectiveAccess(ncFile, ncVarElemGroup);

			int ncVarVrtxSize;
			checkNcError(nc_inq_varid(ncFile, "vertex_size", &ncVarVrtxSize));
			collectiveAccess(ncFile, ncVarVrtxSize);

			int ncVarVrtxCoords;
			checkNcError(nc_inq_varid(ncFile, "vertex_coordinates", &ncVarVrtxCoords));
			collectiveAccess(ncFile, ncVarVrtxCoords);

			Partition* partitions = new Partition[nLocalPart];

			// Read elements
			logInfo(rank) << "Reading elements";
			for (unsigned int i = 0; i < nMaxLocalPart; i++) {
				unsigned int j = i % nLocalPart;

				partitions[j].setRank(j + rank*nMaxLocalPart);

				size_t start[3] = {j + rank*nMaxLocalPart, 0, 0};

				// Element size
				unsigned int size;
				checkNcError(nc_get_var1_uint(ncFile, ncVarElemSize, start, &size));
				partitions[j].setElemSize(size);

				size_t count[3] = {1, size, 4};

				// Elements
				checkNcError(nc_get_vara_int(ncFile, ncVarElemVertices, start, count,
						partitions[j].elements()));
				checkNcError(nc_get_vara_int(ncFile, ncVarElemNeighborRanks, start, count,
						partitions[j].neighborRanks()));
				checkNcError(nc_get_vara_int(ncFile, ncVarElemSideOrientations, start, count,
						partitions[j].neighborOrientation()));
				checkNcError(nc_get_vara_int(ncFile, ncVarElemMPIIndices, start, count,
						partitions[j].mpiIndices()));

				// Boundaries and group
				checkNcError(nc_get_vara_int(ncFile, ncVarElemBoundaries, start, count,
						partitions[j].boundaries()));
				checkNcError(nc_get_vara_int(ncFile, ncVarElemGroup, start, count,
						partitions[j].groups()));

				// Vertex size
				checkNcError(nc_get_var1_uint(ncFile, ncVarVrtxSize, start, &size));
				partitions[j].setVrtxSize(size);

				partitions[j].computeLocalVertices();
			}

			for (unsigned int i = 0; i < nLocalPart; i++) {
				nLocalElements += partitions[i].nElements();
				nLocalVertices += partitions[i].nLocalVertices();
			}

			// Propagate number of local vertices
			unsigned int vertexStart = nLocalVertices;
#ifdef PARALLEL
			MPI_Scan(MPI_IN_PLACE, &vertexStart, 1, MPI_UNSIGNED, MPI_SUM, commIO);
#endif // PARALLEL
			vertexStart -= nLocalVertices;

			logInfo(rank) << "Convert local to global vertex ids";
			for (unsigned int i = 0; i < nLocalPart; i++) {
				partitions[i].convertLocalVertices(vertexStart);
				partitions[i].buildVertexRecvLists();
				vertexStart += partitions[i].nLocalVertices();
			}

			int done;
			do {
				// Transfer vertices until all partitions have received all
				// required global ids. This might take more than one iteration
				// since partitions might have an edge in common but only transfer
				// vertices through faces
				done = true;

				PCU_Comm_Begin();
				for (unsigned int i = 0; i < nLocalPart; i++) {
					std::map<int, std::vector<int>> globalVertexIds;

					partitions[i].buildVertexSendLists(globalVertexIds);

					for (std::map<int, std::vector<int>>::iterator j = globalVertexIds.begin();
							j != globalVertexIds.end(); j++) {

						PCU_Comm_Pack(j->first / nMaxLocalPart, &j->first, sizeof(int));
						int remoteRank = partitions[i].rank();
						PCU_COMM_PACK(j->first / nMaxLocalPart, remoteRank);
						size_t size = j->second.size();
						PCU_COMM_PACK(j->first / nMaxLocalPart, size);
						PCU_Comm_Pack(j->first / nMaxLocalPart, &j->second[0], size*sizeof(int));
					}
				}

				PCU_Comm_Send();

				while (PCU_Comm_Receive()) {
					int localRank;
					PCU_COMM_UNPACK(localRank);
					int remoteRank;
					PCU_COMM_UNPACK(remoteRank);
					size_t size;
					PCU_COMM_UNPACK(size);
					int* globalVertexIds = new int[size];
					PCU_Comm_Unpack(globalVertexIds, size*sizeof(int));

					done &= partitions[localRank % nMaxLocalPart].buildLocal2GlobalMap(remoteRank, globalVertexIds);

					delete [] globalVertexIds;
				}

				for (unsigned int i = 0; i < nLocalPart; i++)
					partitions[i].applyLocal2GlobalMap();

#ifdef PARALLEL
				MPI_Allreduce(MPI_IN_PLACE, &done, 1, MPI_INT, MPI_LAND, commIO);
#endif // PARALLEL
			} while(!done);

			logInfo(rank) << "Collecting elements";
			elements = new int[nLocalElements*4];

			unsigned int pos = 0;
			for (unsigned int i = 0; i < nLocalPart; i++) {
				memcpy(&elements[pos], partitions[i].globalElements(),
						partitions[i].nElements()*4*sizeof(int));
				pos += partitions[i].nElements() * 4;
			}

			logInfo(rank) << "Reading vertices";
			vertices = new double[nLocalVertices*3];

			pos = 0;
			for (unsigned int i = 0; i < nMaxLocalPart; i++) {
				unsigned int j = i % nLocalPart;

				size_t start[3] = {j + rank*nMaxLocalPart, 0, 0};
				size_t count[3] = {1, partitions[j].nVertices(), 3};

				checkNcError(nc_get_vara_double(ncFile, ncVarVrtxCoords, start, count,
						partitions[j].vertices()));

				partitions[j].extractGlobalVertices();

				memcpy(&vertices[pos], partitions[i].globalVertices(),
						partitions[i].nLocalVertices()*3*sizeof(double));
				pos += partitions[i].nLocalVertices() * 3;
			}

			logInfo(rank) << "Collecting boundary conditions";
			boundaries = new int[nLocalElements*4];
			groups = new int[nLocalElements];

			pos = 0;
			for (unsigned int i = 0; i < nLocalPart; i++) {
				partitions[i].convertBoundary();

				memcpy(&boundaries[pos*4], partitions[i].boundaries(),
						partitions[i].nElements()*4*sizeof(int));
				memcpy(&groups[pos], partitions[i].groups(),
						partitions[i].nElements()*sizeof(int));
				pos += partitions[i].nElements();
			}

			delete [] partitions;

			checkNcError(nc_close(ncFile));
		}

		logInfo(rank) << "Constructing the mesh";
		apf::GlobalToVert vertMap;
		apf::construct(m_mesh, elements, nLocalElements, apf::Mesh::TET, vertMap);
		delete [] elements;

		apf::alignMdsRemotes(m_mesh);
		apf::deriveMdsModel(m_mesh);

		logInfo(rank) << "Set coordinates in APF";
		apf::setCoords(m_mesh, vertices, nLocalVertices, vertMap);
		delete [] vertices;

		// Set boundaries
		apf::MeshTag* boundaryTag = m_mesh->createIntTag("boundary condition", 1);
		apf::MeshIterator* it = m_mesh->begin(3);
		unsigned int i = 0;
		while (apf::MeshEntity* element = m_mesh->iterate(it)) {
			apf::Adjacent adjacent;
			m_mesh->getAdjacent(element, 2, adjacent);

			for (unsigned int j = 0; j < 4; j++) {
				if (!boundaries[i*4 + j])
					continue;

				m_mesh->setIntTag(adjacent[j], boundaryTag, &boundaries[i*4 + j]);
			}

			i++;
		}
		m_mesh->end(it);
		delete [] boundaries;

		// Set groups
		apf::MeshTag* groupTag = m_mesh->createIntTag("group", 1);
		it = m_mesh->begin(3);
		i = 0;
		while (apf::MeshEntity* element = m_mesh->iterate(it)) {
			m_mesh->setIntTag(element, groupTag, &groups[i]);
			i++;
		}
		m_mesh->end(it);
		delete [] groups;

#ifdef PARALLEL
		PCU_Switch_Comm(MPI_COMM_WORLD);
#endif // PARALLEL
	}

private:
	/**
	 * Switch to collective access for a netCDf variable
	 */
	static void collectiveAccess(int ncFile, int ncVar)
	{
#ifdef PARALLEL
		checkNcError(nc_var_par_access(ncFile, ncVar, NC_COLLECTIVE));
#endif // PARALLEL
	}

	static void checkNcError(int error)
	{
		if (error != NC_NOERR)
			logError() << "Error while reading netCDF file:" << nc_strerror(error);
	}
};

#endif // NETCDF_MESH_H
