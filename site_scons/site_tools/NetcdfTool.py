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

netcdf_prog_src_serial = """
#include <netcdf.h>

int main() {
    int ncFile;
    nc_open("", 0, &ncFile);
    
    return 0;
}
"""

netcdf_prog_src_parallel = """
#include <netcdf_par.h>

int main() {
    int ncFile;
    nc_open_par("", 0, MPI_COMM_WORLD, MPI_INFO_NULL, &ncFile);
    
    return 0;
}
"""

def CheckNetcdfLinking(context, parallel, message):
    context.Message(message+"... ")
    if parallel:
        src = netcdf_prog_src_parallel
    else:
        src = netcdf_prog_src_serial
    ret = context.TryLink(src, '.c')
    context.Result(ret)
    
    return ret

def generate(env, **kw):
    conf = env.Configure(custom_tests = {'CheckNetcdfLinking' : CheckNetcdfLinking})
    
    if 'parallel' in kw and kw['parallel']:
        parallel = True
        header = 'netcdf_par.h'
    else:
        parallel = False
        header = 'netcdf.h'
        
    if 'required' in kw:
        required = kw['required']
    else:
        required = False
        
    if not conf.CheckLibWithHeader('netcdf', header, 'c'):
        if required:
            print 'Could not find netCDF!'
            env.Exit(1)
        else:
            conf.Finish()
            return
        
    ret = conf.CheckNetcdfLinking(parallel, "Checking whether shared netCDF library is used")
    if not ret:
        # Static library, link with HDF5 and zlib as well
        ret = [conf.CheckLib(lib) for lib in ['hdf5_hl', 'hdf5', 'z']]
        
        if not all(ret):
            if required:
                print 'Could not find HDF5 or zlib!'
                env.Exit(1)
            else:
                conf.Finish()
                return

        # Try to find all additional libraries
        conf.CheckLib('curl')
        conf.CheckLib('gpfs')

        ret = conf.CheckNetcdfLinking(parallel, "Checking whether all netCDF dependencies are found")
        if not ret:
            if required:
                print 'Could not find all netCDF dependencies!'
                env.Exit(1)
            else:
                conf.Finish()
                return
            
    conf.Finish()

def exists(env):
    return True
