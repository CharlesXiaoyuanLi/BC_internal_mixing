#!/bin/python2

### This ncfile is set up to compare the Mass Extinction Coefficients with
### internal vs. external mixing

from netCDF4 import Dataset
import numpy as np
#import scipy.ndimage
#import types
import sys as sys

### Read data files

# Choose type of soot optical property
bctyp = sys.argv[1]
print bctyp

# Open netCDF file
ncfile_im = 'sul-'+bctyp+'-im-ndist.nc'
ncfile_em = 'sul-'+bctyp+'-em-0.0118.nc'
ncfile_cs = 'sul-'+bctyp+'-coat-ndist.new.nc'

fh_im = Dataset(ncfile_im, 'r')
fh_em = Dataset(ncfile_em, 'r')
fh_cs = Dataset(ncfile_cs, 'r')

# Read in coordinates
wl_orig  = fh_im.variables['wl'][:]
vol_orig = fh_im.variables['vol'][:]
volc_orig = fh_cs.variables['vol'][:]
rh_orig  = fh_im.variables['RH'][:]

volpkd = 90
volpkdind = np.where( vol_orig == volpkd )[0].item()
volcpkdind = np.where( volc_orig == volpkd )[0].item()
print volpkdind, volcpkdind
# Read in data
mac_im_orig = fh_im.variables['beta_e'][:][:][:] * (1. -
        np.array(fh_im.variables['ssa'][:][:][:])) 
mac_em_orig = fh_em.variables['beta_e'][:][:][:] * (1. -
        np.array(fh_em.variables['ssa'][:][:][:]))
mac_cs_orig = fh_cs.variables['beta_e'][:][:][:] * (1. -
        np.array(fh_cs.variables['ssa'][:][:][:]))

# Calculate the index
rhpkd = 70
rhpkdind = np.where( rh_orig == rhpkd )[0].item()

# Assign picked values
wlulim = 2
wlulimind = np.where( (wl_orig >= wlulim - 0.1) & (wl_orig <= wlulim + 0.1)
        )[0].item()
wl  = wl_orig[:wlulimind-1:-1]
vol = volpkd
print 'vol = '+str(vol)
rh  = rhpkd
print 'rh = '+str(rh)

mac_im_t = mac_im_orig[:wlulimind-1:-1,volpkdind,:]
mac_im = mac_im_t[:,rhpkdind]
mac_em_t = mac_em_orig[:wlulimind-1:-1,volpkdind,:]
mac_em = mac_em_t[:,rhpkdind]
mac_cs_t = mac_cs_orig[:wlulimind-1:-1,volpkdind,:]
mac_cs = mac_cs_t[:,rhpkdind]
print 'shape = '+str(mac_im.shape)

units = fh_im.variables['beta_e'].units

fh_im.close()
fh_em.close()
fh_cs.close()

# Plot

import matplotlib.pyplot as plt


plt.plot(wl, mac_em, color='r', linewidth=2., label="External")
plt.plot(wl, mac_im, color='b', linewidth=2., label="Homogeneous")
plt.plot(wl, mac_cs, color='g', linewidth=2., label="Core-shell")

plt.ylim([0, 4])
plt.xlim([0, 2.2])
plt.legend()
plt.xlabel('Wavelength ($\mu m$)', fontsize=15)
plt.ylabel('Mass Absorption Coefficient ($m^2/g$)', fontsize=15)

plt.title("BC + Sulfate  ($V_{Sul}/V_{BC+Sul}$ = "+str(volpkd)+"%)")

plt.grid(True)

plt.savefig('mac_wl_'+str(volpkd)+'%_'+bctyp+'.ps',format='ps')
plt.show()
