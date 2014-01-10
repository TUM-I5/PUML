/**
 * @file
 *  This file is part of PUML
 *
 *  For conditions of distribution and use, please see the copyright
 *  notice in the file 'COPYING' at the root directory of this package
 *  and the copyright notice at https://github.com/TUM-I5/PUML
 *
 * @copyright 2013 Technische Universitaet Muenchen
 * @author Sebastian Rettenberger <rettenbs@in.tum.de>
 */

#ifndef GAMBIT_READER_H
#define GAMBIT_READER_H

#include <fstream>
#include <sstream>
#include <vector>

#include "utils/stringutils.h"
#include "utils/logger.h"

/**
 * Describes a boundary face
 */
struct BoundaryFace
{
	/** The element the face belongs to */
	unsigned int element;
	/** The face of the element */
	unsigned int face;
	/** The type of the boundary */
	unsigned int type;
};

class GambitReader
{
private:
	struct GambitSection {
		/** Start of the section */
		size_t seekPosition;
		/** Number of elements (lines) in the section */
		unsigned int nLines;
		/** Line size */
		size_t lineSize;
	};

	struct ElementSection : GambitSection {
		/** Start of vertices in an element */
		size_t vertexStart;
		/** Size per vertex id */
		size_t vertexSize;
	};

	struct GroupSection : GambitSection {
		/** Size per element id */
		size_t elementSize;
		/** Material type */
		int material;
	};

	struct BoundarySection : GambitSection {
		/** Type of the boundary */
		int type;
	};

	std::ifstream m_mesh;

	// Section descriptions
	GambitSection m_vertices;
	ElementSection m_elements;
	std::vector<GroupSection> m_groups;
	std::vector<BoundarySection> m_boundaries;

public:
	GambitReader()
	{
	}

	GambitReader(const char* meshFile)
	{
		open(meshFile);
	}

	void open(const char* meshFile)
	{
		m_mesh.open(meshFile);

		// Files are readable?
		if (!m_mesh)
			logError() << "Could not open mesh file" << meshFile;

		std::string line;

		// Read header information
		getline(m_mesh, line);	// First line contains version
								// we ignore this line for know
		line.clear();
		getline(m_mesh, line);
		utils::StringUtils::trim(line);
		if (line != GAMBIT_FILE_ID)
			logError() << "Not a Gambit mesh file:" << meshFile;

		std::string name;
		getline(m_mesh, name);		// Internal name

		getline(m_mesh, line); 	// Skip line:
								// PROGRAM: Gambit VERSION: x.y.z

		getline(m_mesh, line);		// Date
		//trim(line);

		getline(m_mesh, line); 	// Skip problem size names

		unsigned int nGroups, nBoundaries, dimensions;
		m_mesh >> m_vertices.nLines;
		m_mesh >> m_elements.nLines;
		m_mesh >> nGroups;
		m_mesh >> nBoundaries;
		m_mesh >> dimensions;
		getline(m_mesh, line); 	// Skip rest of the line
		if (dimensions != 3)
			logError() << "Gambit file does not contain a 3 dimensional mesh";

		getline(m_mesh, line);
		utils::StringUtils::trim(line);
		if (line != ENDSECTION)
			logError() << "Invalid Gambit format: End of header expected, found" << line;

		// Find seek positions
		// Vertices
		getline(m_mesh, line);
		if (line.find(NODAL_COORDINATES) == std::string::npos)
			logError() << "Invalid Gambit format: Coordinates expected, found" << line;
		m_vertices.seekPosition = m_mesh.tellg();
		getline(m_mesh, line);
		m_vertices.lineSize = line.length() + 1;
		m_mesh.seekg(m_vertices.seekPosition + m_vertices.nLines * m_vertices.lineSize);
		getline(m_mesh, line);
		utils::StringUtils::rtrim(line); // remove \r
		if (line != ENDSECTION)
			logError() << "Invalid Gambit format: End of coordinates expected, found" << line;

		// Elements
		getline(m_mesh, line);
		if (line.find(ELEMENT_CELLS) == std::string::npos)
			logError() << "Invalid Gambit format: Elements expected, found" << line;
		m_elements.seekPosition = m_mesh.tellg();
		getline(m_mesh, line);
		m_elements.lineSize = line.length() + 1;
		m_mesh.seekg(m_elements.seekPosition + m_elements.nLines * m_elements.lineSize);

		std::istringstream ss(line);
		int tmp;
		ss >> tmp; // Element id
		ss >> tmp;
		ss >> tmp;
		ss.seekg(1, std::stringstream::cur); // White space
		m_elements.vertexStart = ss.tellg();
		// TODO works only for tetrahedrals
		utils::StringUtils::rtrim(line);
		if ((line.length() - m_elements.vertexStart) % 4 != 0)
			logError() << "Invalid Gambit format: Mesh does not seem to contain tetrahedrals";
		m_elements.vertexSize = (line.length() - m_elements.vertexStart) / 4;

		getline(m_mesh, line);
		utils::StringUtils::rtrim(line); // remove \r
		if (line != ENDSECTION)
			logError() << "Invalid Gambit format: End of elements expected, found" << line;

		// Groups
		m_groups.resize(nGroups);
		for (unsigned int i = 0; i < nGroups; i++) {
			getline(m_mesh, line);
			if (line.find(ELEMENT_GROUP) == std::string::npos)
				logError() << "Invalid Gambit format: Group expected, found" << line;
			getline(m_mesh, line);

			// Group size and material
			std::string y;
			std::istringstream ss(line);
			ss >> y;
			ss >> m_groups[i].material;
			ss >> y;
			ss >> m_groups[i].nLines;	// This is not the actual number of lines
										// Because Gambit stores more than one element
										// per line.

			getline(m_mesh, line);
			getline(m_mesh, line);
			m_groups[i].seekPosition = m_mesh.tellg();

			getline(m_mesh, line);
			m_groups[i].lineSize = line.size() + 1;

			utils::StringUtils::rtrim(line);
			if (line.length() % ELEMENTS_PER_LINE_GROUP != 0)
				logError() << "Invalid Gambit format: Mesh does not contain" << ELEMENTS_PER_LINE_GROUP << "in one group line";
			m_groups[i].elementSize = line.length() / ELEMENTS_PER_LINE_GROUP;

			m_mesh.seekg(m_groups[i].seekPosition
					+ (m_groups[i].nLines / ELEMENTS_PER_LINE_GROUP) * m_groups[i].lineSize);

			if (m_groups[i].nLines % ELEMENTS_PER_LINE_GROUP != 0)
				// Last line
				m_mesh.seekg(m_groups[i].lineSize - (ELEMENTS_PER_LINE_GROUP
						- m_groups[i].nLines % ELEMENTS_PER_LINE_GROUP) * m_groups[i].elementSize,
						std::ifstream::cur);

			getline(m_mesh, line);
			utils::StringUtils::rtrim(line); // remove \r
			if (line != ENDSECTION)
				logError() << "Invalid Gambit format: End of group expected, found" << line;
		}

		// Boundaries
		m_boundaries.resize(nBoundaries);
		for (unsigned int i = 0; i < nBoundaries; i++) {
			getline(m_mesh, line);
			if (line.find(BOUNDARY_CONDITIONS) == std::string::npos)
				logError() << "Invalid Gambit format: Boundaries expected, found" << line;
			getline(m_mesh, line);
			m_boundaries[i].seekPosition = m_mesh.tellg();

			// Boundary type and size
			int x;
			std::istringstream ss(line);
			ss >> m_boundaries[i].type;
			ss >> x;
			ss >> m_boundaries[i].nLines;

			getline(m_mesh, line);
			m_boundaries[i].lineSize = line.size() + 1;

			m_mesh.seekg(m_boundaries[i].seekPosition
					+ m_boundaries[i].nLines * m_boundaries[i].lineSize);

			getline(m_mesh, line);
			utils::StringUtils::rtrim(line); // remove \r
			if (line != ENDSECTION)
				logError() << "Invalid Gambit format: End of boundaries expected, found" << line;
		}
	}

	unsigned int nVertices()
	{
		return m_vertices.nLines;
	}

	unsigned int nElements()
	{
		return m_elements.nLines;
	}

	/**
	 * @return Number of boundary faces
	 */
	unsigned int nBoundaries()
	{
		unsigned int count = 0;
		for (std::vector<BoundarySection>::const_iterator i = m_boundaries.begin();
				i != m_boundaries.end(); i++) {
			count += i->nLines;
		}

		return count;
	}

	/**
	 * Reads vertices from start tp start+count from the file and stores them in vertices
	 *
	 * @param vertices Buffer to store coordinates of each vertex. The caller is responsible
	 *  for allocating the buffer. The size of the buffer must be count*dimensions
	 *
	 * @todo Only 3 dimensional meshes are supported
	 */
	void readVertices(unsigned int start, unsigned int count, double* vertices)
	{
		m_mesh.clear();

		m_mesh.seekg(m_vertices.seekPosition + start * m_vertices.lineSize
				+ m_vertices.lineSize - 3*COORDINATE_SIZE - 1);

		char* buf = new char[3*COORDINATE_SIZE];

		for (unsigned int i = 0; i < count; i++) {
			m_mesh.read(buf, 3*COORDINATE_SIZE);

			for (int j = 0; j < 3; j++) {
				std::istringstream ss(std::string(&buf[j*COORDINATE_SIZE], COORDINATE_SIZE));
				ss >> vertices[i * 3 + j];
			}

			m_mesh.seekg(m_vertices.lineSize - 3*COORDINATE_SIZE, std::fstream::cur);
		}

		delete [] buf;
	}

	/**
	 * Reads elements from start to start+count from the file and stores them in elements
	 *
	 * @param elements Buffer to store vertices of each element. The caller is responsible
	 *  for allocating the buffer. The Size of the buffer must be count*vertices_per_element.
	 *
	 * @todo Only tetrahedral meshes are supported
	 * @todo Support for varying coordinate/vertexid fields
	 */
	void readElements(unsigned int start, unsigned int count, unsigned int* elements)
	{
		m_mesh.clear();

		m_mesh.seekg(m_elements.seekPosition + start * m_elements.lineSize + m_elements.vertexStart);

		char* buf = new char[4*m_elements.vertexSize];

		for (unsigned int i = 0; i < count; i++) {
			m_mesh.read(buf, 4*m_elements.vertexSize);

			for (int j = 0; j < 4; j++) {
				std::istringstream ss(std::string(&buf[j*m_elements.vertexSize], m_elements.vertexSize));
				ss >> elements[i * 4 + j];
				elements[i * 4 + j]--;
			}

			m_mesh.seekg(m_elements.lineSize - 4*m_elements.vertexSize, std::fstream::cur);
		}

		delete [] buf;
	}

	/**
	 * Reads boundaries from start to start+count from the file and stores them in
	 * <code>boundaries</code>.
	 *
	 * @param elements Buffer to store boundaries. Each boundary consists of the element, the
	 *  face and the boundary type.
	 */
	void readBoundaries(unsigned int start, unsigned int count, BoundaryFace* boundaries)
	{
		m_mesh.clear();

		// Get the boundary, were we should start reading
		std::vector<BoundarySection>::const_iterator section;
		for (section = m_boundaries.begin();
				section != m_boundaries.end() && section->nLines < start; section++)
			start -= section->nLines;

		m_mesh.seekg(section->seekPosition + start * section->lineSize);

		char* buf = new char[ELEMENT_SIZE_BOUNDARY + ELEMENT_TYPE + FACE_SIZE];

		for (unsigned int i = 0; i < count; i++) {
			m_mesh.read(buf, ELEMENT_SIZE_BOUNDARY + ELEMENT_TYPE + FACE_SIZE);

			std::istringstream ssE(std::string(buf, ELEMENT_SIZE_BOUNDARY));
			ssE >> boundaries[i].element;
			boundaries[i].element--;

			std::istringstream ssF(std::string(&buf[ELEMENT_SIZE_BOUNDARY + ELEMENT_TYPE], FACE_SIZE));
			ssF >> boundaries[i].face;
			boundaries[i].face = face2internal(boundaries[i].face);

			boundaries[i].type = section->type;

			start++; // Line in the current section
			if (start >= section->nLines) {
				start = 0;
				section++;

				m_mesh.seekg(section->seekPosition);
			} else {
				m_mesh.seekg(section->lineSize - (ELEMENT_SIZE_BOUNDARY + ELEMENT_TYPE + FACE_SIZE),
						std::fstream::cur);
			}
		}

		delete [] buf;
	}

private:
	/**
	 * Convert the gambit face number to the internal face number
	 */
	static unsigned int face2internal(unsigned int face)
	{
		switch (face) {
		case 1:
		case 2:
			return face-1;
		case 3:
			return 3;
		case 4:
			return 2;
		}

		return -1;
	}

	/** Number of character required to store a coordinate */
	static const size_t COORDINATE_SIZE = 20ul;
	/** Number of elements stored in one group line */
	static const size_t ELEMENTS_PER_LINE_GROUP = 10ul;
	/** Number of characters required to store an element id in boundary conditions */
	static const size_t ELEMENT_SIZE_BOUNDARY = 10ul;
	/** Number of characters required to store an element type */
	static const size_t ELEMENT_TYPE = 5ul;
	/** Number of characters required to store a face id */
	static const size_t FACE_SIZE = 5ul;

	static const char* GAMBIT_FILE_ID;
	static const char* ENDSECTION;
	static const char* NODAL_COORDINATES;
	static const char* ELEMENT_CELLS;
	static const char* ELEMENT_GROUP;
	static const char* BOUNDARY_CONDITIONS;
};

#endif // GAMBIT_READER_H
