#! /usr/bin/python

# @file
# This file is part of SeisSol.
# 
# @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
# 
# 
# @section LICENSE
# This software was developed at Technische Universitaet Muenchen, who is the owner of the software.
# 
# According to good scientific practice, publications on results achieved in whole or in part due to this software should cite at least one paper or referring to an URL presenting the this software software.
# 
# The owner wishes to make the software available to all users to use, reproduce, modify, distribute and redistribute also for commercial purposes under the following conditions of the original BSD license. Linking this software module statically or dynamically with other modules is making a combined work based on this software. Thus, the terms and conditions of this license cover the whole combination. As a special exception, the copyright holders of this software give you permission to link it with independent modules or to instantiate templates and macros from this software's source files to produce an executable, regardless of the license terms of these independent modules, and to copy and distribute the resulting executable under terms of your choice, provided that you also meet, for each linked independent module, the terms and conditions of this license of that module.
# 
# Copyright (c) 2013
# Technische Universitaet Muenchen
# Department of Informatics
# Chair of Scientific Computing
# http://www5.in.tum.de/
# 
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
# 
# Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# All advertising materials mentioning features or use of this software must display the following acknowledgement: This product includes software developed by the Technische Universitaet Muenchen (TUM), Germany, and its contributors.
# Neither the name of the Technische Universitaet Muenchen, Munich, Germany nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
#
# @section DESCRIPTION
#
# Finds netCDF and add inlcude pathes, libs and library pathes
#

netcdf_prog_src = """
#include HEADER

int main() {
    int ncFile;
    nc_open("", 0, &ncFile);
    
    return 0;
}
"""

def CheckNetcdfLinking(context, header, message):
    context.Message(message+"... ")
    ret = context.TryLink(netcdf_prog_src.replace('HEADER', '<'+header+'>', 1), '.c')
    context.Result(ret)
    
    return ret

def generate(env, **kw):
    conf = env.Configure(custom_tests = {'CheckNetcdfLinking' : CheckNetcdfLinking})
    
    if 'parallel' in kw and kw['parallel']:
        header = 'netcdf_par.h'
    else:
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
            return
        
    ret = conf.CheckNetcdfLinking(header, "Checking whether shared netCDF library is used")
    if not ret:
        # Static library, link with HDF5 and zlib as well
        ret = [conf.CheckLib(lib) for lib in ['hdf5_hl', 'hdf5', 'z']]
        
        if not all(ret):
            if required:
                print 'Could not find HDF5 or zlib!'
                env.Exit(1)
            else:
                return

        # Try to find all additional libraries
        conf.CheckLib('curl')
        conf.CheckLib('gpfs')

        ret = conf.CheckNetcdfLinking(header, "Checking wheater all netCDF dependencies are found")
        if not ret:
            if required:
                print 'Could not find all netCDF dependencies!'
                env.Exit(1)
            else:
                return

def exists(env):
    return True
