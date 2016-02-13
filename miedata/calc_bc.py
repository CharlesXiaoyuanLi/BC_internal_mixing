#!/bin/python2

### This ncfile is set up to compare the Mass Extinction Coefficients with
### internal vs. external mixing

from netCDF4 import Dataset
import numpy as np
import sys as sys
import itertools

for bctype in ("","m","m_lim1","h","hh"):

    print 'bc'+bctype
### Read data files

# Open netCDF file
    ncfile_bc = 'bc'+bctype+'.nc'

    fh_bc = Dataset(ncfile_bc, 'r')

# Read in coordinates
    wl_orig  = fh_bc.variables['wl'][:]
    rm_orig = fh_bc.variables['RadMean'][:]

    wlpkd = 0.55  # choose the wavelength at 550nm to represent shortwave radiation
    wlpkdind = np.where((wl_orig >= wlpkd-0.02) & (wl_orig <= wlpkd+0.02))[0].item()
#    print wlpkdind

# Read in data
    mac_bc_orig = fh_bc.variables['beta_e'][:][:] * (1. -
            np.array(fh_bc.variables['ssa'][:][:])) 
    ssa_bc_orig = fh_bc.variables['ssa'][:][:]

# Assign picked values
    wl  = wlpkd

    mac_bc = mac_bc_orig[wlpkdind,:]
    ssa_bc = ssa_bc_orig[wlpkdind,:]
#    print 'shape = '+str(mac_bc.shape)

    units = fh_bc.variables['beta_e'].units

    for i in range(len(rm_orig)):
        print 'Mean Radius = '+str(rm_orig[i])
        print 'Pure BC MAC = '+str(mac_bc[i])+' '+units
        print 'Pure BC SSA = '+str(ssa_bc[i])


    fh_bc.close()
