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

#define private public

#include "PUML/NetcdfEntity.h"
#include "PUML/NetcdfGroup.h"
#include "PUML/NetcdfPum.h"

static const char* TEST_FILENAME = "test.nc.pum";

class TestNetcdfEntity : public CxxTest::TestSuite
{
private:
	PUML::NetcdfPum m_ncPum;
	PUML::NetcdfGroup* m_ncGroup;
	PUML::NetcdfEntity* m_ncEntity0;
	PUML::NetcdfEntity* m_ncEntity1;

public:
	void setUp()
	{
#ifdef PARALLEL
		TS_ASSERT(m_ncPum.create(TEST_FILENAME, 5, MPI_COMM_WORLD));
#else // PARALLEL
		TS_ASSERT(m_ncPum.create(TEST_FILENAME, 2));
#endif // PARALLEL
		//m_ncGroup = m_ncPum.createGroup("testGroup", 25);
		m_ncGroup = m_ncPum.createGroup("testGroup");
		TS_ASSERT(m_ncGroup);

		// Without extra dimension
		m_ncEntity0 = m_ncGroup->createEntity("testEntity0", PUML::Type::Int);
		TS_ASSERT(m_ncEntity0);

		// With extra dimension
		PUML::Dimension dim = m_ncGroup->createDimension("testDimension", 2);
		m_ncEntity1 = m_ncGroup->createEntity("testEntity1", PUML::Type::Float, 1, &dim);
		TS_ASSERT(m_ncEntity1);

		TS_ASSERT(m_ncPum.endDefinition());

		int r = 0;
		int s = 1;
#ifdef PARALLEL
		MPI_Comm_rank(MPI_COMM_WORLD, &r);
		MPI_Comm_size(MPI_COMM_WORLD, &s);
#endif // PARALLEL

		// 2 MPI processes / 5 partitions
		TS_ASSERT(m_ncGroup->setSize(r, 5));
		TS_ASSERT(m_ncGroup->setSize(r+s, 5));
		TS_ASSERT(m_ncGroup->setSize(r+2*s, 5));
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

	void testSetCollective()
	{
		TS_ASSERT(m_ncEntity0->setCollective(true));
	}

	void testPut()
	{
		TS_ASSERT(m_ncEntity0->setCollective(true));
		TS_ASSERT(m_ncEntity1->setCollective(true));

		int r = 0;
#ifdef PARALLEL
		MPI_Comm_rank(MPI_COMM_WORLD, &r);
#endif // PARALLEL

		float values[2*5];
		for (int i = 0; i < 2*5; i++)
			values[i] = i+1000*r;

		TS_ASSERT(m_ncEntity0->put(r, 5, values));
		TS_ASSERT(m_ncEntity1->put(r, 5, values));
	}

	void testGet()
	{
		TS_ASSERT(m_ncEntity0->setCollective(true));
		TS_ASSERT(m_ncEntity1->setCollective(true));

		int r = 0;
#ifdef PARALLEL
		MPI_Comm_rank(MPI_COMM_WORLD, &r);
#endif // PARALLEL

		float values[2*5];
		for (int i = 0; i < 2*5; i++)
			values[i] = i+1000*r;

		TS_ASSERT(m_ncEntity0->put(r, 5, values));
		TS_ASSERT(m_ncEntity1->put(r, 5, values));

		setUpOpen();

		TS_ASSERT(m_ncEntity0->setCollective(true));
		TS_ASSERT(m_ncEntity1->setCollective(true));

		for (int i = 0; i < 2*5; i++)
			values[i] = 0;
		TS_ASSERT(m_ncEntity1->get(r, 5, values));
		for (int i = r; i < r+5; i++)
			TS_ASSERT_EQUALS(values[i], i+1000*r);

		for (int i = 0; i < 2*5; i++)
			values[i] = 0;
		TS_ASSERT(m_ncEntity1->get(r, values));
		for (int i = r; i < r+5; i++)
			TS_ASSERT_EQUALS(values[i], i+1000*r);
		// TODO check why we can't get ints as floats
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

		m_ncEntity0 = m_ncGroup->getEntity("testEntity0");
		m_ncEntity1 = m_ncGroup->getEntity("testEntity1");
	}
};
