#!/bin/python2

### This ncfile is set up to compare the Mass Extinction Coefficients with
### internal vs. external mixing of three species: BC, sulfate and OC.

from netCDF4 import Dataset
import numpy as np
#import scipy.ndimage
#import types
import sys as sys

### Read data files

# Open netCDF file
ncfile_im = 'bcm-sul-oc-coat.nc'
ncfile_em = 'bcm-sul-coat-oc-em.nc'

fh_im = Dataset(ncfile_im, 'r')
fh_em = Dataset(ncfile_em, 'r')

# Read in coordinates
wl_orig  = fh_im.variables['wl'][:]
vols_orig = fh_im.variables['vols'][:]
volo_orig = fh_im.variables['volo'][:]
rh_orig  = fh_im.variables['RH'][:]

wlpkd = 0.55  # choose the wavelength at 550nm to represent shortwave radiation
wlpkdind = np.where((wl_orig >= wlpkd-0.02) & (wl_orig <= wlpkd+0.02))[0].item()

# Read in data
mec_im_orig = fh_im.variables['beta_e'][:][:][:][:]
mec_em_orig = fh_em.variables['beta_e'][:][:][:][:]

ssa_im_orig = fh_im.variables['ssa'][:][:][:][:]
ssa_em_orig = fh_em.variables['ssa'][:][:][:][:]

g_im_orig = fh_im.variables['g'][:][:][:][:]
g_em_orig = fh_em.variables['g'][:][:][:][:]


# Calculate the index
rhpkd = np.array([30,40,50,60,70,80,84,88,91])
rhpkdind = [np.where( rh_orig == s )[0].item() for s in rhpkd]

volspkd = 90
volspkdind = np.where( vols_orig == volspkd )[0].item()

volopkd = np.array([70,75,80,84,88,92,96,98])
volopkdind = [np.where( volo_orig == s )[0].item() for s in volopkd]

print 'vol_OC = '+str(volo_orig[volopkdind])
print 'rh = '+str(rh_orig[rhpkdind])

# Assign picked values
mec_im_t = mec_im_orig[wlpkdind,volspkdind,volopkdind,:]
mec_im = mec_im_t[:,rhpkdind]
mec_em_t = mec_em_orig[wlpkdind,volspkdind,volopkdind,:]
mec_em = mec_em_t[:,rhpkdind]

units = fh_im.variables['beta_e'].units

ssa_im_t = ssa_im_orig[wlpkdind,volspkdind,volopkdind,:]
ssa_im = ssa_im_t[:,rhpkdind]
ssa_em_t = ssa_em_orig[wlpkdind,volspkdind,volopkdind,:]
ssa_em = ssa_em_t[:,rhpkdind]

g_im_t = g_im_orig[wlpkdind,volspkdind,volopkdind,:]
g_im = g_im_t[:,rhpkdind]
g_em_t = g_em_orig[wlpkdind,volspkdind,volopkdind,:]
g_em = g_em_t[:,rhpkdind]


fh_im.close()
fh_em.close()

# Plot

import matplotlib.pyplot as plt

vol = volo_orig[volopkdind]
rh  = rh_orig[rhpkdind]

area_t = 15 * (vol/45.)**2
area = np.zeros((vol.size,rh.size),float)
for i in np.arange(rh.size):
    area[:,i] = area_t

color_t = rh
colors = np.zeros((vol.size,rh.size),float)
for i in np.arange(vol.size):
    colors[i,:] = color_t

fig = plt.figure()

#fig,axes = plt.subplots(nrows=3,ncols=2)

ax1=plt.subplot(3,1,1)
ax1.scatter(mec_em, mec_im, c=colors, cmap='RdYlBu', s=area, alpha=0.8,
        linewidths=0.5, marker='o')

#ax1.set_xlim([2, 18])
#ax1.set_ylim([2, 18])

lims = [
    np.min([plt.xlim(), plt.ylim()]),
    np.max([plt.ylim(), plt.ylim()]),
    ]

ax1.plot(lims,lims,'k-',alpha=0.75,zorder=0)

ax1.set_ylabel('Homogeneous IM', fontsize=8)
ax1.set_title('Mass Extinction Coefficient ($m^2/g$)', fontsize=10)
ax1.set_xlabel('External Mixing', fontsize=8)
ax1.set_aspect('equal')

plt.grid(True)
plt.tick_params(labelsize=7)


plt.subplot(3,1,2)
plt.scatter(ssa_em, ssa_im, c=colors, cmap='RdYlBu', s=area, alpha=0.8,
        linewidths=0.5, marker='o')
#plt.xlim([.3, 1.])
#plt.ylim([.3, 1.])
lims = [
    np.min([plt.xlim(), plt.ylim()]),
    np.max([plt.ylim(), plt.ylim()]),
    ]
plt.plot(lims,lims,'k-',alpha=0.75,zorder=0)
plt.ylabel('Homogeneous IM', fontsize=8)
plt.xlabel('External Mixing', fontsize=8)
plt.title('Single Scattering Albedo', fontsize=10)
plt.grid(True)
plt.tick_params(labelsize=7)
ax3 = plt.subplot(3,1,2)
ax3.set_aspect('equal')


plt.subplot(3,1,3)
im=plt.scatter(g_em, g_im, c=colors, cmap='RdYlBu', s=area, alpha=0.8,
        linewidths=0.5, marker='o')
#plt.xlim([0.55, 0.8])
#plt.ylim([0.55, 0.8])
lims = [
    np.min([plt.xlim(), plt.ylim()]),
    np.max([plt.ylim(), plt.ylim()]),
    ]
plt.plot(lims,lims,'k-',alpha=0.75,zorder=0)
plt.ylabel('Homogeneous IM', fontsize=8)
plt.xlabel('External Mixing', fontsize=8)
plt.title('Asymmetry Factor', fontsize=10)
plt.grid(True)
plt.tick_params(labelsize=7)
ax5 = plt.subplot(3,1,3)
ax5.set_aspect('equal')


plt.tight_layout()

fig.subplots_adjust(bottom=0.13)
cax = fig.add_axes([0.13, 0.05, 0.8, 0.015])
cbar = fig.colorbar(im,cax=cax,ticks=rhpkd,orientation='horizontal')
cbar.solids.set_edgecolor("face")
cbar.set_label("RH (%)")


plt.savefig('radppty_bcm-sul-im-oc.pdf',format='pdf')
plt.show()
