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

#include "FidapReader.h"

const char* FidapReader::FIDAP_FILE_ID = "** FIDAP NEUTRAL FILE";
const char* FidapReader::NODAL_COORDINATES = "NODAL COORDINATES";
const char* FidapReader::BOUNDARY_CONDITIONS = "BOUNDARY CONDITIONS";
const char* FidapReader::ELEMENT_GROUPS = "ELEMENT GROUPS";

const char* FidapReader::ZONE_GROUP = "GROUP:";
const char* FidapReader::ZONE_ELEMENTS = "ELEMENTS:";
const char* FidapReader::ZONE_NODES = "NODES:";
const char* FidapReader::ZONE_GEOMETRY = "GEOMETRY:";
const char* FidapReader::ZONE_TYPE = "TYPE:";
