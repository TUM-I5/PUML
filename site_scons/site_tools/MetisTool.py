#! /usr/bin/python

# @file
#  This file is part of PUML
#
#  For conditions of distribution and use, please see the copyright
#  notice in the file 'COPYING' at the root directory of this package
#  and the copyright notice at https://github.com/TUM-I5/PUML
# 
# @copyright 2013 Technische Universitaet Muenchen
# @author Sebastian Rettenberger <rettenbs@in.tum.de>
#

def generate(env, **kw):
    conf = env.Configure()
        
    if 'parallel' in kw and kw['parallel']:
        header = 'parmetis.h'
        lib = 'parmetis'
    else:
        header = 'metis.h'
        lib = 'metis'
        
    if 'required' in kw:
        required = kw['required']
    else:
        required = False
        
    if not conf.CheckLibWithHeader(lib, header, 'c'):
        if required:
            print 'Could not find METIS!'
            env.Exit(1)
        else:
            conf.Finish()
            return
            
    # For the parallel version we need to link against metis as well
    if 'parallel' in kw and kw['parallel']:
        if not conf.CheckLib('metis'):
            if required:
                print 'Could not find ParMETIS'
                env.Exit(1)
            else:
                conf.Finish()
                return
            
    conf.Finish()

def exists(env):
    return True