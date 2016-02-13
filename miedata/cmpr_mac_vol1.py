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
ncfile_im = 'sul-'+bctyp+'-im.nc'
ncfile_em = 'sul-'+bctyp+'-em.nc'
ncfile_cs = 'sul-'+bctyp+'-coat.nc'

fh_im = Dataset(ncfile_im, 'r')
fh_em = Dataset(ncfile_em, 'r')
fh_cs = Dataset(ncfile_cs, 'r')

# Read in coordinates
wl_orig  = fh_im.variables['wl'][:]
vol_orig = fh_im.variables['vol'][:]
rh_orig  = fh_im.variables['RH'][:]

wlpkd = 0.55  # choose the wavelength at 550nm to represent shortwave radiation
wlpkdind = np.where((wl_orig >= wlpkd-0.02) & (wl_orig <= wlpkd+0.02))[0].item()
print wlpkdind

# Read in data
mac_im_orig = fh_im.variables['beta_e'][:][:][:] * (1. -
        np.array(fh_im.variables['ssa'][:][:][:])) 
mac_em_orig = fh_em.variables['beta_e'][:][:][:] * (1. -
        np.array(fh_em.variables['ssa'][:][:][:]))
mac_cs_orig = fh_cs.variables['beta_e'][:][:][:] * (1. -
        np.array(fh_cs.variables['ssa'][:][:][:]))

# Calculate the index
rhpkd = np.array([30,40,50,60,70,80,84,88,91,93,96])
rhpkdind = [np.where( rh_orig == s )[0].item() for s in rhpkd]

# Assign picked values
wl  = wlpkd
volint = -1
vol = vol_orig[::volint]
print 'vol = '+str(vol)
rh  = rhpkd
print 'rh = '+str(rh)

mac_im_t = mac_im_orig[wlpkdind,::volint,:]
mac_im = mac_im_t[:,rhpkdind]
mac_em_t = mac_em_orig[wlpkdind,::volint,:]
mac_em = mac_em_t[:,rhpkdind]
mac_cs_t = mac_cs_orig[wlpkdind,::volint,:]
mac_cs = mac_cs_t[:,rhpkdind]
print 'shape = '+str(mac_im.shape)

units = fh_im.variables['beta_e'].units

print bctyp.upper()+' Pure BC IM mac = '+str(mac_im[0,0])+' '+units
print vol_orig[-1]
print rh_orig[0]
print bctyp.upper()+' Pure BC SSA = '+str(1.-fh_im.variables['ssa'][wlpkdind][-1][0])

fh_im.close()
fh_em.close()
fh_cs.close()

print mac_im[:,0]

# Plot

import matplotlib.pyplot as plt

plt.plot(100-vol, mac_em[:,10], color='r', linewidth=1, label='EM')
plt.plot(100-vol, mac_im[:,10], color='b', linewidth=1, label='IM')
plt.plot(100-vol, mac_cs[:,10], color='g', linewidth=1, label='CS')

plt.legend()

plt.ylim([0, 12])

plt.xlabel('BC VMR (%)', fontsize=15)
plt.ylabel('Mass Absorption Coefficient ($m^2/g$)', fontsize=15)

plt.grid(True)

#plt.savefig('mac_'+bctyp+'.pdf',format='pdf')
plt.show()
