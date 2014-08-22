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

#ifndef APF_NATIVE_H
#define APF_NATIVE_H

#include <apfMDS.h>
#include <gmi_mesh.h>
#include <gmi_null.h>

#include "utils/logger.h"

#include "MeshInput.h"

class ApfNative : public MeshInput
{
public:
	ApfNative(const char* mesh, const char* model = 0L)
	{
		if (model)
			gmi_register_mesh();
		else {
			gmi_register_null();
			model = ".null";
		}
		m_mesh = apf::loadMdsMesh(model, mesh);
	}
};

#endif // APF_NATIVE_H
