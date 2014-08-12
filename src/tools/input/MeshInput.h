/**
 * @file
 *  This file is part of PUML
 *
 *  For conditions of distribution and use, please see the copyright
 *  notice in the file 'COPYING' at the root directory of this package
 *  and the copyright notice at https://github.com/TUM-I5/PUML
 *
 * @copyright 2014 Technische Universitaet Muenchen
 * @author Sebastian Rettenberger <rettenbs@in.tum.de>
 */

#ifndef MESH_INTPUT_H
#define MESH_INTPUT_H

/**
 * Interface for mesh input
 */
class MeshInput
{
private:
	unsigned int m_nVertices;
	unsigned int m_nElements;

	/** Number of local vertices */
	unsigned int m_nLocalVertices;
	/** Number of local elements */
	unsigned int m_nLocalElements;
public:
	virtual ~MeshInput() {}

	unsigned int nVertices() const
	{
		return m_nVertices;
	}

	unsigned int nElements() const
	{
		return m_nElements;
	}

	/**
	 * Number of vertices this process is responsible for
	 */
	unsigned int nLocalVertices() const
	{
		return m_nLocalVertices;
	}

	unsigned int nLocalElements() const
	{
		return m_nLocalElements;
	}

	virtual int rankOfVert(unsigned int vertex) const = 0;
	virtual unsigned int posOfVert(unsigned int vertex) const = 0;

	virtual unsigned int elemStart(int rank) const = 0;
	virtual int rankOfElem(unsigned int element) const = 0;
	virtual unsigned int posOfElem(unsigned int element) const = 0;

	virtual void getVertices(double* vertices) = 0;
	virtual void getElements(unsigned int* elements) = 0;
	virtual void getGroups(unsigned int* groups) = 0;
	virtual void getBoundaries(unsigned int* boundaries) = 0;

protected:
	void setNVertices(unsigned int nVertices)
	{
		m_nVertices = nVertices;
	}

	void setNElements(unsigned int nElements)
	{
		m_nElements = nElements;
	}

	void setNLocalVertices(unsigned int nLocalVertices)
	{
		m_nLocalVertices = nLocalVertices;
	}

	void setNLocalElements(unsigned int nLocalElements)
	{
		m_nLocalElements = nLocalElements;
	}
};

#endif // MESH_INPUT_H
