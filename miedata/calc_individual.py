#!/bin/python2

### This ncfile is set up to compare the Mass Extinction Coefficients with
### internal vs. external mixing

from netCDF4 import Dataset
import numpy as np
import sys as sys
import itertools

### Read data files

# Open netCDF file
ncfile_bc = 'bc_size.nc'
ncfile_su = 'sulfate_size.nc'

fh_bc = Dataset(ncfile_bc, 'r')
fh_su = Dataset(ncfile_su, 'r')

# Read in coordinates
wl_orig  = fh_bc.variables['wl'][:]
rm_orig = fh_bc.variables['RadMean'][:]
rh_orig  = fh_bc.variables['RH'][:]

wlpkd = 0.55  # choose the wavelength at 550nm to represent shortwave radiation
wlpkdind = np.where((wl_orig >= wlpkd-0.02) & (wl_orig <= wlpkd+0.02))[0].item()
print wlpkdind

# Read in data
mac_bc_orig = fh_bc.variables['beta_e'][:][:][:] * (1. -
        np.array(fh_bc.variables['ssa'][:][:][:])) 
ssa_bc_orig = fh_bc.variables['ssa'][:][:][:]
mec_su_orig = fh_su.variables['beta_e'][:][:][:]
ssa_su_orig = fh_su.variables['ssa'][:][:][:]

# Assign picked values
wl  = wlpkd

mac_bc = mac_bc_orig[wlpkdind,:,:]
ssa_bc = ssa_bc_orig[wlpkdind,:,:]
mec_su = mec_su_orig[wlpkdind,:,:]
ssa_su = ssa_su_orig[wlpkdind,:,:]
print 'shape = '+str(mac_bc.shape)

units = fh_bc.variables['beta_e'].units

for i in range(len(rm_orig)):
    print 'Mean Radius = '+str(rm_orig[i])
    print 'Pure BC MAC = '+str(mac_bc[i][0])+' '+units
    print 'Pure BC SSA = '+str(ssa_bc[i][0])
for i,j in itertools.product(range(len(rm_orig)),range(len(rh_orig))):
    print 'RH = '+str(rh_orig[j])+'%'
    print 'Pure Sulfate MEC = '+str(mec_su[i][j])+' '+units
    print 'Pure Sulfate SSA = '+str(ssa_su[i][j])


fh_bc.close()
fh_su.close()
