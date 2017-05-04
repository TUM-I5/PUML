#! /usr/bin/python

# @file
#  This file is part of PUML
#
#  For conditions of distribution and use, please see the copyright
#  notice in the file 'COPYING' at the root directory of this package
#  and the copyright notice at https://github.com/TUM-I5/PUML
#
# @copyright 2014-2016 Technische Universitaet Muenchen
# @author Sebastian Rettenberger <rettenbs@in.tum.de>
#

import os

def tryLibPath(env, conf, libPath, mpiWrapper, setRpath):
    oldLibPath = env['LIBPATH']
    oldRPath = env['RPATH']

    success = False

    try:
        env.AppendUnique(LIBPATH=libPath)

        # Add path for parasolid library
        psLibPath = [p for p in map(lambda p: os.path.join(p, 'psKrnl'), env['LIBPATH']) if os.path.exists(p)]
        env.AppendUnique(LIBPATH=psLibPath)
        if setRpath:
            env.AppendUnique(RPATH=psLibPath)

        # TODO check all headers
        # TODO not all libraries may be available/required
        libs = [
            ('SimAdvMeshing', 'SimAdvMeshing.h'),
            ('SimMeshing', 'MeshSim.h'),
            ('SimField', 'SimField.h'),
            ('SimExport', 'SimExport.h'),
            ('SimDiscrete', 'SimDiscrete.h'),
            ('SimMeshTools', 'SimMeshTools.h'),
            (['SimParasolid260', 'SimParasolid270', 'SimParasolid280'], None),
            ('SimPartitionedMesh-mpi', 'SimPartitionedMesh.h'),
            ('SimPartitionWrapper-'+mpiWrapper, None),
            ('SimModel', 'SimModel.h'),
            ('pskernel', None)
        ]

        for l in libs:
            if l[1]:
                if not conf.CheckLibWithHeader(l[0], l[1], 'c++'):
                    return False
            else:
                if not conf.CheckLib(l[0]):
                    return False

        success = True
    finally:
        if not success:
            env['LIBPATH'] = oldLibPath
            env['RPATH'] = oldRPath

    return True


def generate(env, **kw):
    conf = env.Configure()

    if 'required' in kw:
        required = kw['required']
    else:
        required = False

    if 'rpath' in kw:
        rpath = kw['rpath']
    else:
        rpath = True

    if 'mpiLib' in kw:
        mpi = kw['mpiLib']
    else:
        mpi = 'mpich2'

    # TODO support different configurations
    for lib in ['x64_rhel5_gcc41', 'x64_rhel6_gcc44']:
        libPath = [p for p in map(lambda p: os.path.join(p, lib), env['LIBPATH']) if os.path.exists(p)]
        if not libPath:
            # This library directory does not exist -> try next
            continue

        print 'Checking for simmodeler libraries in ' + lib + '...'
        if tryLibPath(env, conf, libPath, mpi, rpath):
            conf.Finish()
            return

    # Not found
    if required:
        print 'Could not find SimModSuite!'
        env.Exit(1)

    conf.Finish()

def exists(env):
    return True
