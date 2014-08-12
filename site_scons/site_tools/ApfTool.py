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

def generate(env, **kw):
    conf = env.Configure()
        
    if 'required' in kw:
        required = kw['required']
    else:
        required = False
        
    if not conf.CheckLibWithHeader('gmi_sim', 'apf.h', 'c++'):
        if required:
            print 'Could not find APF!'
            env.Exit(1)
        else:
            conf.Finish()
            return
        
    # Additional libs
    libs = ['gmi', 'mds', 'apf_sim', 'apf', 'pcu']
    for l in libs:
        if not conf.CheckLib(l):
            if required:
                print 'Could not find APF. ' + l + ' is missing!'
                env.Exit(1)
            
    conf.Finish()

def exists(env):
    return True