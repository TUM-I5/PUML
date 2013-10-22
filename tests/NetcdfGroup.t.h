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

#ifdef PARALLEL
#include <mpi.h>
#endif // PARALLEL

#include <cstdio>

#include <cxxtest/TestSuite.h>

#include "PUML/NetcdfEntity.h"
#include "PUML/NetcdfGroup.h"
#include "PUML/NetcdfPum.h"

static const char* TEST_FILENAME = "test.nc.pum";

class TestNetcdfGroup : public CxxTest::TestSuite
{
private:
	PUML::NetcdfPum m_ncPum;
	PUML::NetcdfGroup* m_ncGroup;
	PUML::NetcdfGroup* m_ncIndexedGroup;

public:
	void setUp()
	{
#ifdef PARALLEL
		TS_ASSERT(m_ncPum.create(TEST_FILENAME, 5, MPI_COMM_WORLD));
#else // PARALLEL
		TS_ASSERT(m_ncPum.create(TEST_FILENAME, 2));
#endif // PARALLEL
		m_ncGroup = m_ncPum.createGroup("testGroup");
		TS_ASSERT(m_ncGroup);

		m_ncIndexedGroup = m_ncPum.createGroupIndexed("testIndexedGroup");
		TS_ASSERT(m_ncIndexedGroup);
	}

	void tearDown()
	{
		if (!m_ncPum.isValid())
			TS_FAIL(m_ncPum.errorMsg());
		TS_ASSERT(m_ncPum.close());

		// Remove generated test file
		remove(TEST_FILENAME);
	}

	/**
	 * Test the load constructor only
	 */
	void testConstructor()
	{
		setUpOpen();
	}

	void testCreateDimension()
	{
		m_ncGroup->createDimension("testDim", 2);
		TS_ASSERT(m_ncGroup->isValid());
	}

	void testCreateEntity()
	{
		PUML::NetcdfEntity* e = m_ncGroup->createEntity("testEnt", PUML::Type::Int);
		TS_ASSERT(e);

		// Entity with another dimension
		PUML::Dimension d = m_ncGroup->createDimension("testDim", 2);
		TS_ASSERT(m_ncGroup->isValid());
		e = m_ncGroup->createEntity("testEntWithDim", PUML::Type::Int, 1, &d);
		TS_ASSERT(e);
	}

	void testSetSize()
	{
		TS_ASSERT(m_ncPum.endDefinition());

		int r = 0;
		int s = 1;
#ifdef PARALLEL
		MPI_Comm_rank(MPI_COMM_WORLD, &r);
		MPI_Comm_size(MPI_COMM_WORLD, &s);
#endif // PARALLEL

		// Not the first partition -> fail
		TS_ASSERT(!m_ncGroup->setSize(r+1, 5));

		TS_ASSERT(m_ncGroup->setSize(r, 5));
		TS_ASSERT(m_ncGroup->setSize(r+s, 5));
		TS_ASSERT(m_ncGroup->setSize(r+2*s, 5));
	}

	void testSize()
	{
		testSetSize();

		setUpOpen();

		TS_ASSERT_EQUALS(m_ncGroup->size(0), 5ul);
#ifdef PARALLEL
		// _size is not set because we did not insert any values
		TS_ASSERT_EQUALS(m_ncGroup->size(1), 5ul);
#endif // PARALLEL
	}

	void testSetIndex()
	{

		TS_ASSERT(m_ncPum.endDefinition());

		int r = 0;
		int s = 1;
#ifdef PARALLEL
		MPI_Comm_rank(MPI_COMM_WORLD, &r);
		MPI_Comm_size(MPI_COMM_WORLD, &s);
#endif // PARALLEL

		// Not the first partition -> fail
		TS_ASSERT(!m_ncIndexedGroup->setSize(r+1, 5));

		TS_ASSERT(m_ncIndexedGroup->setSize(r, 5));
		TS_ASSERT(m_ncIndexedGroup->setSize(r+s, 5));
		TS_ASSERT(m_ncIndexedGroup->setSize(r+2*s, 5));

		unsigned long index[] = {r, r+1, r+2, 2, 3};
		TS_ASSERT(m_ncIndexedGroup->putIndex(r, 5, index));
	}

private:
	void setUpOpen()
	{
		TS_ASSERT(m_ncPum.close());

#ifdef PARALLEL
		TS_ASSERT(m_ncPum.open(TEST_FILENAME, MPI_COMM_WORLD));
#else // PARALLEL
		TS_ASSERT(m_ncPum.open(TEST_FILENAME));
#endif // PARALLEL

		m_ncGroup = m_ncPum.getGroup("testGroup");
	}
};
