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

import os
import sys

#
# set possible variables
#
vars = Variables()

# read parameters from a file if given
vars.AddVariables(
  PathVariable( 'buildVariablesFile', 'location of the python file, which contains the build variables', None, PathVariable.PathIsFile )
)
env = Environment(variables=vars)
if 'buildVariablesFile' in env:
  vars = Variables(env['buildVariablesFile'])

# PUML specific variables
vars.AddVariables(
  PathVariable( 'buildDir', 'where to build the code', 'build', PathVariable.PathIsDirCreate ),

  EnumVariable( 'compileMode', 'mode of the compilation', 'release',
                allowed_values=('debug', 'release')
              ),

  EnumVariable( 'parallelization', 'level of parallelization', 'none',
                allowed_values=('none', 'mpi')
              ),

  BoolVariable( 'unitTests', 'builds additional unit tests',
                False
              ),

  EnumVariable( 'logLevel',
                'logging level. \'debug\' prints all information available, \'info\' prints information at runtime (time step, plot number), \'warning\' prints warnings during runtime, \'error\' is most basic and prints errors only',
                'info',
                allowed_values=('debug', 'info', 'warning', 'error')
              ),

  BoolVariable( 'simModSuite', 'compile with support for simModSuite from Simmetrix',
                False
              ),

  ( 'mpiLib', 'MPI library against this is linked (only required for simModSuite)',
                'mpich2',
                None, None
              )
)

env.Tool('PrefixPathTool')

# external variables
vars.AddVariables(
  env['PREFIX_PATH_VARIABLE'],

  PathVariable( 'cc',
                'C compiler (default: gcc (serial), mpicc (parallel))',
                None,
                PathVariable.PathAccept ),

  PathVariable( 'cxx',
                'C++ compiler (default: g++ (serial), mpiCC (parallel))',
                None,
                PathVariable.PathAccept ),

  BoolVariable( 'useEnv',
                'set variables set in the execution environment',
                True )
)

# generate help text
Help(vars.GenerateHelpText(env))
if '-h' in sys.argv or '--help' in sys.argv:
  import SCons
  print SCons.Script.help_text
  env.Exit()

# handle unknown, maybe misspelled variables
unknownVariables = vars.UnknownVariables()

# remove the buildVariablesFile from the list of unknown variables (used before)
if 'buildVariablesFile' in unknownVariables:
  unknownVariables.pop('buildVariablesFile')

# exit in the case of unknown variables
if unknownVariables:
  print "*** The following build variables are unknown: " + str(unknownVariables.keys())
  env.Exit(1)

# set environment
env = Environment(variables=vars)

# Set environment
if env['useEnv']:
  env['ENV'] = os.environ

#
# precompiler, compiler and linker flags
#

# set compiler
if 'cc' in env:
  env['CC'] = env['cxx']
else:
  if env['parallelization'] in ['mpi']:
    env['CC'] = 'mpicc'
  else:
    env['CC'] = 'gcc'

if 'cxx' in env:
  env['CXX'] = env['cxx']
else:
  if env['parallelization'] in ['mpi']:
    env['CXX'] = 'mpicxx'
  else:
    env['CXX'] = 'g++'

# add parallel flag for mpi
if env['parallelization'] in ['mpi']:
  env.Append(CPPDEFINES=['PARALLEL'])

# set level of logger
if env['logLevel'] == 'debug':
  env.Append(CPPDEFINES=['LOG_LEVEL=3'])
elif env['logLevel'] == 'info':
  env.Append(CPPDEFINES=['LOG_LEVEL=2'])
elif env['logLevel'] == 'warning':
  env.Append(CPPDEFINES=['LOG_LEVEL=1'])
elif env['logLevel'] == 'error':
  env.Append(CPPDEFINES=['LOG_LEVEL=0'])
else:
  assert(false)

# compiler flags for generated kernels
env.Append(CXXFLAGS = ['-Wall', '-ansi', '-std=c++0x'])
if env['compileMode'] == 'debug':
  env.Append(CXXFLAGS=['-O0','-g'])
elif env['compileMode'] == 'release':
  env.Append(CPPDEFINES=['NDEBUG'])
  env.Append(CXXFLAGS=['-O2'])

# add pathname to the list of directories which are search for include
env.Append(CPPPATH=['#/src'])

# Add prefix path
env.Tool('PrefixPathTool')

# netCDF
env.Tool('NetcdfTool', parallel=(env['parallelization'] in ['mpi']), required=True)

# proj.4
env.Tool('ProjTool', required=False)

#
# setup the library name and the build directory
#
lib_name = 'puml'

# Parallelization?
if env['parallelization'] != 'mpi':
    lib_name += '_nompi'

# build directory
env['libFile'] = env['buildDir']+'/'+lib_name
env['execDir'] = env['buildDir']+'/bin'
env['buildDir'] = env['buildDir']+'/build_'+lib_name

# setup additional dependencies for tools
env.tools = env.Clone()
env.tools.Append(CPPPATH=['#/submodules'])
# TODO does not work (at the moment it is not required)
#env.tools.Append(LIBS=[os.path.split(env['libFile'])[1]])
#env.tools.Append(LIBPATH=[os.path.split(env['libFile'])[0]])
env.tools.Append(CXXFLAGS = ['-fopenmp'])
env.tools.Append(LINKFLAGS= ['-fopenmp'])
# (Par)METIS
env.tools.Tool('MetisTool', parallel=(env['parallelization'] in ['mpi']), required=True)

if env['parallelization'] in ['mpi']:
    # APF
    env.tools.Tool('ApfTool', required=True, sim=env['simModSuite'])

    # SimModSuite
    if env['simModSuite']:
        env.tools.Tool('SimModSuiteTool', mpiLib=env['mpiLib'], required=True)
        env.tools.Append(CPPDEFINES=['USE_SIMMOD'])

# get the source files
env.sourceFiles = {}
env.sourceFiles['lib'] = []

Export('env')
SConscript('src/SConscript', variant_dir='#/'+env['buildDir'], src_dir='#/', duplicate=0)
Import('env')

# build standard version
env.StaticLibrary('#/'+env['libFile'], env.sourceFiles['lib'])

# build tools
if env['parallelization'] in ['mpi']:
    env.tools.Program('#/'+env['execDir']+'/pumgen', env.sourceFiles['pumgen'])

# build unit tests
if env['unitTests']:
    # Use our own env for testing
    env = env.Clone()

    # Some MPI version requires this flag
    if env['parallelization'] == 'mpi':
        env.Append(CPPDEFINES=['MPICH_IGNORE_CXX_SEEK'])

    # define location of cxxtest
    env['CXXTEST'] = 'submodules/cxxtest'

    # Continue testing if tests fail
    env['CXXTEST_SKIP_ERRORS'] = True

    # Parallel tests?
    if env['parallelization'] in ['mpi']:
        env['CXXTEST_COMMAND'] = 'mpiexec -np 2 %t'
        env['CXXTEST_OPTS'] = ['--main=cxxtest_main']

    # add cxxtest-tool
    env.Tool('cxxtest', toolpath=[env['CXXTEST']+'/build_tools/SCons'])

    # Get test source files
    env.sourceFiles['tests'] = []

    Export('env')
    SConscript('tests/SConscript', variant_dir='#/'+env['buildDir']+'/tests', src_dir='#/')
    Import('env')

    # build unit tests
    env.CxxTest(target='#/'+env['buildDir']+'/tests/cxxtest_runner',
                source=env.sourceFiles['lib']+env.sourceFiles['tests'])
