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

#include "GambitReader.h"

const char* GambitReader::GAMBIT_FILE_ID = "** GAMBIT NEUTRAL FILE";
const char* GambitReader::ENDSECTION = "ENDOFSECTION";
const char* GambitReader::NODAL_COORDINATES = "NODAL COORDINATES";
const char* GambitReader::ELEMENT_CELLS = "ELEMENTS/CELLS";
const char* GambitReader::ELEMENT_GROUP = "ELEMENT GROUP";
const char* GambitReader::BOUNDARY_CONDITIONS = "BOUNDARY CONDITIONS";
