#!/bin/python2

### This ncfile is set up to compare the Mass Extinction Coefficients with
### internal vs. external mixing

from netCDF4 import Dataset
import numpy as np
#import scipy.ndimage
#import types
import sys as sys

### Read data files

# Open netCDF file
ncfile_em = 'bc_size.nc'

fh_em = Dataset(ncfile_em, 'r')

# Read in coordinates
wl_orig  = fh_em.variables['wl'][:]
rad_orig = fh_em.variables['RadMean'][:]
rh_orig  = fh_em.variables['RH'][:]

# Read in data
mac_em_orig = fh_em.variables['beta_e'][:][:][:] * (1. -
        np.array(fh_em.variables['ssa'][:][:][:]))

# Calculate the index
rhpkdind = 0

# Assign picked values
wlulim = 3
wlulimind = np.where( (wl_orig >= wlulim - 0.1) & (wl_orig <= wlulim + 0.1)
        )[0].item()
wl  = wl_orig[:wlulimind-1:-1]

mac_em_t = mac_em_orig[:wlulimind-1:-1,:,:]
print 'shape = '+str(mac_em_t.shape)
mac_em = mac_em_t[:,:,rhpkdind]
print 'shape = '+str(mac_em.shape)

units = fh_em.variables['beta_e'].units

fh_em.close()

# Plot

import matplotlib.pyplot as plt


plt.plot(wl, mac_em[:,0], color='r', linewidth=2., label="$D_{pg} = 0.0118 \mu m$")
plt.plot(wl, mac_em[:,1], color='b', linewidth=2., label="$D_{pg} = 0.0355 \mu m$")
plt.plot(wl, mac_em[:,2], color='g', linewidth=2., label="$D_{pg} = 0.05 \mu m$")

#plt.ylim([0, 5])
plt.legend()
plt.xlabel('Wavelength ($\mu m$)', fontsize=15)
plt.ylabel('Mass Absorption Coefficient ($m^2/g$)', fontsize=15)

plt.grid(True)

plt.savefig('mac_bc_wl.pdf',format='pdf')
plt.show()
