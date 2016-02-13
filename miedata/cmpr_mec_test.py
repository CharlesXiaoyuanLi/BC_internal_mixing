### This ncfile is set up to compare the Mass Extinction Coefficients with
### internal vs. external mixing

from netCDF4 import Dataset
import numpy as np
#import scipy.ndimage
#import types
import sys as sys

### Read data files

# Open netCDF file
ncfile_im = 'sul-bc-im.nc'
ncfile_em = 'sul-bc-em.nc'
ncfile_cs = 'sul-bc-coat.nc'

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
mec_im_orig = fh_im.variables['beta_e'][:][:][:]
mec_em_orig = fh_em.variables['beta_e'][:][:][:]
mec_cs_orig = fh_cs.variables['beta_e'][:][:][:]

# Calculate the index
rhpkd = np.array([50,60,70,80,90])
rhpkdind = [np.where( rh_orig == s )[0].item() for s in rhpkd]
print rhpkdind

# Assign picked values
wl  = wlpkd
vol = vol_orig
rh  = rhpkd

mec_im_t = mec_im_orig[wlpkdind,:,:]
mec_im = mec_im_t[:,rhpkdind]
mec_em_t = mec_em_orig[wlpkdind,:,:]
mec_em = mec_em_t[:,rhpkdind]
mec_cs_t = mec_cs_orig[wlpkdind,:,:]
mec_cs = mec_cs_t[:,rhpkdind]
print mec_im.shape

units = fh_im.variables['beta_e'].units

fh_im.close()
fh_em.close()
fh_cs.close()

# Plot

import matplotlib.pyplot as plt

#area_t = 15 * ((1.-vol)/20)**2
area_t = 30
area = np.zeros((vol.size,rh.size),float)
for i in np.arange(rh.size):
    area[:,i] = area_t

color_t = rh
colors = np.zeros((vol.size,rh.size),float)
for i in np.arange(vol.size):
    colors[i,:] = color_t

fig, ax = plt.subplots()

ax.scatter(mec_em, mec_cs, c=colors, cmap='Spectral', s=area, alpha=0.1,
        linewidths=0.5, marker='^')

#ax.set_xscale('log')
#ax.set_yscale('log')

ax.set_xlim([4, 14])
ax.set_ylim([4, 14])

ax.set_xlabel('External Mixing', fontsize=15)
ax.set_ylabel('Homogeneous Internal Mixing', fontsize=15)
ax.set_title('Mass Extinction Coefficient (BC + Sulfate)')

ax.grid(True)
fig.tight_layout()

plt.savefig('mec_test.ps',format='ps')
