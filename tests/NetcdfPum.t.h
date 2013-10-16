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

#include <cxxtest/GlobalFixture.h>
#include <cxxtest/TestSuite.h>

#include "PUML/NetcdfGroup.h"
#include "PUML/NetcdfPum.h"

#ifdef PARALLEL
int cxxtest_main(int, char**);

static bool mpiInitSuccess = true;

int main(int argc, char** argv)
{
	if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
		mpiInitSuccess = false;

	int exitCode = cxxtest_main(argc, argv);

	MPI_Finalize();

	return exitCode;
}

/**
 * Init and finalize MPI (requires the main functions)
 *
 * @see main
 */
class MPIFixture : public CxxTest::GlobalFixture
{
public:
	bool setUpWorld()
	{
		return mpiInitSuccess;
	}
};
static MPIFixture mpiFixture;
#endif // PARALLEL

static const char* TEST_FILENAME = "test.nc.pum";

class TestNetcdfPum : public CxxTest::TestSuite
{
private:
	PUML::NetcdfPum m_ncPum;

public:
	void setUp()
	{
#ifdef PARALLEL
		TS_ASSERT(m_ncPum.create(TEST_FILENAME, 2, MPI_COMM_WORLD));
#else // PARALLEL
		TS_ASSERT(m_ncPum.create(TEST_FILENAME, 2));
#endif // PARALLEL
	}

	void tearDown()
	{
		if (!m_ncPum.isValid())
			TS_FAIL(m_ncPum.errorMsg());
		TS_ASSERT(m_ncPum.close());

		// Remove generated test file
		remove(TEST_FILENAME);
	}

	void testCreate()
	{
		TS_ASSERT(m_ncPum.isValid());
	}

	void testOpen()
	{
		setUpOpen();
#ifdef PARALLEL
		TS_ASSERT(m_ncPum.open(TEST_FILENAME, MPI_COMM_WORLD));
#else // PARALLEL
		TS_ASSERT(m_ncPum.open(TEST_FILENAME));
#endif // PARALLEL
	}

	void testCreateGroup()
	{
		PUML::NetcdfGroup ncGroup = m_ncPum.createGroup("testGroup");
		TS_ASSERT(ncGroup.isValid());
	}

	void testEndDefinition()
	{
		PUML::NetcdfGroup ncGroup = m_ncPum.createGroup("testGroup");
		TS_ASSERT(m_ncPum.endDefinition());
	}

private:
	void setUpOpen()
	{
		TS_ASSERT(m_ncPum.close());
	}
};
