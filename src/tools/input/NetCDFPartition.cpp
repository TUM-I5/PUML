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

#include "NetCDFPartition.h"

const int Partition::FACE2VERTICES[4][3] = {
	{0, 2, 1},
	{0, 1, 3},
	{0, 3, 2},
	{1, 2, 3}
};

const int Partition::INTERNAL2EX_ORDER[4] = {0, 1, 3, 2};
