#! /usr/bin/python

# @file
#  This file is part of PUML
#
#  For conditions of distribution and use, please see the copyright
#  notice in the file 'COPYING' at the root directory of this package
#  and the copyright notice at https://github.com/TUM-I5/PUML
# 
# @copyright 2014 Technische Universitaet Muenchen
# @author Sebastian Rettenberger <rettenbs@in.tum.de>
#

import os

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
    libPath = [p for p in map(lambda p: os.path.join(p, 'x64_rhel5_gcc41'), env['LIBPATH']) if os.path.exists(p)]
    env.AppendUnique(LIBPATH=libPath)
    
    # TODO check which headers belong to which library
    if not conf.CheckLibWithHeader('SimMeshing', 'MeshSim.h', 'c++'):
        if required:
            print 'Could not find SimModSuite!'
            env.Exit(1)
        else:
            conf.Finish()
            return
    
    # Add path for parasolid library    
    libPath = [p for p in map(lambda p: os.path.join(p, 'psKrnl'), env['LIBPATH']) if os.path.exists(p)]
    env.AppendUnique(LIBPATH=libPath)
    if rpath:
        env.AppendUnique(RPATH=libPath)
    
    # TODO not all libraries may be available/required
    # TODO different parasolid versions not handled currently
    libs = ['SimMeshTools', 'SimParasolid260', 'SimPartitionedMesh-mpi',
            'SimPartitionWrapper-'+mpi, 'SimModel', 'pskernel']
    for l in libs:
        if not conf.CheckLib(l):
            if required:
                print 'Could not find SimModSuite. ' + l + ' is missing!'
                env.Exit(1)
    
    conf.Finish()
    
def exists(env):
    return True