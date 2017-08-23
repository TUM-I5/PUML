

from netCDF4 import Dataset
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Transformed Asagi file into a netcdf readable by PUMGEN (--velocity-model AsagiDerived)')
parser.add_argument('AsagiFile', help='Asagi input filename')
args = parser.parse_args()


#read ASAGI file
fname = args.AsagiFile
fh = Dataset(fname, 'r', format='NETCDF3_64BIT')
x = fh.variables['x'][:]
y = fh.variables['y'][:]
z = fh.variables['z'][:]

nx = x.shape[0]
ny = y.shape[0]
nz = z.shape[0]

encodedkeys = fh.variables.keys()
keys = [str(encodedkeys[x]) for x in range(len(encodedkeys))]
if 'data' in keys:
   data = fh.variables['data'][:]
   #compute Vp
   Vp = np.zeros(data['rho'].shape)
   Vp = np.sqrt((data['lambda'] + 2.0 * data['mu'])/data['rho'])
elif all(x in keys for x in ['rho', 'g', 'lambda']):
   rho = fh.variables['rho'][:,:,:]
   mu = fh.variables['g'][:,:,:]
   lambd = fh.variables['lambda'][:,:,:]
   Vp = np.zeros(rho.shape)
   Vp = np.sqrt((lambd + 2.0 * mu)/rho)
else:
   print "unknown netcdf structure..."
   print "data or (rho,g,lambda) not in netcdf file"
   print keys
   exit()   

fh.close()
print "min max Vp: %f %f" %(np.amin(Vp), np.amax(Vp))
#write Asagi derived file

fname = 'AsagiDerivedVp.nc'
fnetcdf = Dataset(fname, 'w', format='NETCDF3_64BIT')
fnetcdf.createDimension('x', nx)
fnetcdf.createDimension('y', ny)
fnetcdf.createDimension('z', nz)
fnetcdf.createDimension('xyz', nx*ny*nz)
desc = 'New Zealand 3d velocity model: only Vp'
fnetcdf.description = desc

x_nc = fnetcdf.createVariable('x', 'd', ('x'))
y_nc = fnetcdf.createVariable('y', 'd', ('y'))
z_nc = fnetcdf.createVariable('z', 'd', ('z'))

Vp_nc = fnetcdf.createVariable('Vp', 'd', ('xyz'))

x_nc[:] = x
y_nc[:] = y
z_nc[:] = z

Vp_nc[:] = Vp.flatten()

fnetcdf.close()


