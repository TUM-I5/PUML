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
        
    if 'sim' in kw:
        sim = kw['sim']
    else:
        sim = False
        
    # Libs required to couple SimModSuite
    simLibs = [('gmi_sim', 'gmi_sim.h'), ('apf_sim', 'apfSIM.h')]
    # Other APF libs
    libs = [('gmi', 'gmi.h'), ('mds', 'apfMDS.h'), ('ma', 'ma.h'),
            ('apf', 'apf.h'), ('pcu', 'PCU.h')]
    
    if sim:
        libs = simLibs + libs
        
    for lib in libs:
        if not conf.CheckLibWithHeader(lib[0], lib[1], 'c++'):
            if required:
                print 'Could not find APF!'
                env.Exit(1)
            
    conf.Finish()

def exists(env):
    return True