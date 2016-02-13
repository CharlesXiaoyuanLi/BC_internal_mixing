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
ncfile_im = './bcm-sul-oc-im.nc'
ncfile_em = './bcm-sul-oc-em.nc'
ncfile_cs = './bcm-sul-oc-coat.nc'

fh_im = Dataset(ncfile_im, 'r')
fh_em = Dataset(ncfile_em, 'r')
fh_cs = Dataset(ncfile_cs, 'r')

# Read in coordinates
wl_orig  = fh_im.variables['wl'][:]
vols_orig = fh_im.variables['vols'][:]
volo_orig = fh_im.variables['volo'][:]
rh_orig  = fh_im.variables['RH'][:]

volspkd = 90
volspkdind = np.where( vols_orig == volspkd )[0].item()
volopkd = 90
volopkdind = np.where( volo_orig == volopkd )[0].item()

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
wlulimind = np.where( (wl_orig >= wlulim - 0.1) & (wl_orig <= wlulim + 0.1) )[0].item()
wl  = wl_orig[:wlulimind-1:-1]
print 'wl = '+str(wl)
vols = volspkd
print 'vols = '+str(vols)
volo = volopkd
print 'volo = '+str(volo)
rh  = rhpkd
print 'rh = '+str(rh)

mac_im_t = mac_im_orig[:wlulimind-1:-1,volspkdind,volopkdind,:]
mac_im = mac_im_t[:,rhpkdind]
mac_em_t = mac_em_orig[:wlulimind-1:-1,volspkdind,volopkdind,:]
mac_em = mac_em_t[:,rhpkdind]
mac_cs_t = mac_cs_orig[:wlulimind-1:-1,volspkdind,volopkdind,:]
mac_cs = mac_cs_t[:,rhpkdind]
print 'shape = '+str(mac_im.shape)

units = fh_im.variables['beta_e'].units

fh_im.close()
fh_em.close()
fh_cs.close()

print mac_em
print mac_im
print mac_cs


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

plt.title("BC + Sulfate + OC  ($V_{Sul}/V_{BC+Sul}$ = "+str(volspkd)+"%, $V_{OC}/V_{BC+OC}$ = "+str(volopkd)+"%)")

plt.grid(True)

plt.savefig('mac_wl_'+str(volspkd)+'%s_'+str(volopkd)+'%o.ps',format='ps')
plt.show()
