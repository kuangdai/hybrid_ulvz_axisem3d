import numpy as np
from netCDF4 import Dataset

# Define the filename
filename = 'ulvz.nc'

# Create a new NetCDF file
with Dataset(filename, 'w', format='NETCDF4') as ncfile:
    # Create dimensions
    ncfile.createDimension('radius', 2)
    ncfile.createDimension('theta', 2)
    ncfile.createDimension('phi', 2)

    # Create variables
    radius = ncfile.createVariable('radius', 'f4', ('radius',))
    theta = ncfile.createVariable('theta', 'f4', ('theta',))
    phi = ncfile.createVariable('phi', 'f4', ('phi',))
    vs = ncfile.createVariable('vs', 'f4', ('radius', 'theta', 'phi'))
    vp = ncfile.createVariable('vp', 'f4', ('radius', 'theta', 'phi'))
    rho = ncfile.createVariable('rho', 'f4', ('radius', 'theta', 'phi'))

    # Assign data to variables
    radius[:] = [3480, 3510]
    theta[:] = [0, 8.57142857143]
    phi[:] = [0, 360]

    vs[:, :, :] = -0.2
    vp[:, :, :] = -0.1
    rho[:, :, :] = 0.05

# Print message confirming file creation
print(f"NetCDF file '{filename}' created successfully.")