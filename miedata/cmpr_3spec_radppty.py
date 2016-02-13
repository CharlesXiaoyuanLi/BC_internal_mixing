#!/bin/python2

### This ncfile is set up to compare the Mass Extinction Coefficients with
### internal vs. external mixing

from netCDF4 import Dataset
import numpy as np
#import scipy.ndimage
#import types
import sys as sys

### Read data files

volspkd = 10

# Open netCDF file
ncfile_3im = 'bcm-sul-oc-im.nc'
ncfile_2im = 'bcm-sul-im-oc-em.nc'
ncfile_3cs = 'bcm-sul-oc-coat.nc'
ncfile_2cs = 'bcm-sul-coat-oc-em.nc'

fh_3im = Dataset(ncfile_3im, 'r')
fh_2im = Dataset(ncfile_2im, 'r')
fh_3cs = Dataset(ncfile_3cs, 'r')
fh_2cs = Dataset(ncfile_2cs, 'r')

# Read in coordinates
wl_orig  = fh_3im.variables['wl'][:]
vols_orig = fh_3im.variables['vols'][:]
volo_orig = fh_3im.variables['volo'][:]
rh_orig  = fh_3im.variables['RH'][:]

wlpkd = 0.55  # choose the wavelength at 550nm to represent shortwave radiation
wlpkdind = np.where((wl_orig >= wlpkd-0.02) & (wl_orig <= wlpkd+0.02))[0].item()

# Read in data
mec_3im_orig = fh_3im.variables['beta_e'][:][:][:]
mec_2im_orig = fh_2im.variables['beta_e'][:][:][:]
mec_3cs_orig = fh_3cs.variables['beta_e'][:][:][:]
mec_2cs_orig = fh_2cs.variables['beta_e'][:][:][:]

ssa_3im_orig = fh_3im.variables['ssa'][:][:][:]
ssa_2im_orig = fh_2im.variables['ssa'][:][:][:]
ssa_3cs_orig = fh_3cs.variables['ssa'][:][:][:]
ssa_2cs_orig = fh_2cs.variables['ssa'][:][:][:]

g_3im_orig = fh_3im.variables['g'][:][:][:]
g_2im_orig = fh_2im.variables['g'][:][:][:]
g_3cs_orig = fh_3cs.variables['g'][:][:][:]
g_2cs_orig = fh_2cs.variables['g'][:][:][:]


# Calculate the index
rhpkd = np.array([30,60,80,90,94])
rhpkdind = [np.where( rh_orig == s )[0].item() for s in rhpkd]

volspkdind = np.where( vols_orig == volspkd )[0].item()

volopkd = np.array([70,80,90,94,98])
volopkdind = [np.where( volo_orig == s )[0].item() for s in volopkd]

print 'vol_OC = '+str(volo_orig[volopkdind])
print 'vol_sulfate = '+str(vols_orig[volspkdind])
print 'rh = '+str(rh_orig[rhpkdind])

mec_3im_t = mec_3im_orig[wlpkdind,volspkdind,volopkdind,:]
mec_3im = mec_3im_t[:,rhpkdind]
mec_2im_t = mec_2im_orig[wlpkdind,volspkdind,volopkdind,:]
mec_2im = mec_2im_t[:,rhpkdind]
mec_3cs_t = mec_3cs_orig[wlpkdind,volspkdind,volopkdind,:]
mec_3cs = mec_3cs_t[:,rhpkdind]
mec_2cs_t = mec_2cs_orig[wlpkdind,volspkdind,volopkdind,:]
mec_2cs = mec_2cs_t[:,rhpkdind]

units = fh_3im.variables['beta_e'].units

ssa_3im_t = ssa_3im_orig[wlpkdind,volspkdind,volopkdind,:]
ssa_3im = ssa_3im_t[:,rhpkdind]
ssa_2im_t = ssa_2im_orig[wlpkdind,volspkdind,volopkdind,:]
ssa_2im = ssa_2im_t[:,rhpkdind]
ssa_3cs_t = ssa_3cs_orig[wlpkdind,volspkdind,volopkdind,:]
ssa_3cs = ssa_3cs_t[:,rhpkdind]
ssa_2cs_t = ssa_2cs_orig[wlpkdind,volspkdind,volopkdind,:]
ssa_2cs = ssa_2cs_t[:,rhpkdind]

g_3im_t = g_3im_orig[wlpkdind,volopkdind,volspkdind,:]
g_3im = g_3im_t[:,rhpkdind]
g_2im_t = g_2im_orig[wlpkdind,volopkdind,volspkdind,:]
g_2im = g_2im_t[:,rhpkdind]
g_3cs_t = g_3cs_orig[wlpkdind,volopkdind,volspkdind,:]
g_3cs = g_3cs_t[:,rhpkdind]
g_2cs_t = g_2cs_orig[wlpkdind,volopkdind,volspkdind,:]
g_2cs = g_2cs_t[:,rhpkdind]


fh_3im.close()
fh_2im.close()
fh_3cs.close()
fh_2cs.close()

# Plot

import matplotlib.pyplot as plt

vol = volo_orig[volopkdind]
rh  = rh_orig[rhpkdind]

area_t = 15 * (vol/40.)**2
area = np.zeros((vol.size,rh.size),float)
for i in np.arange(rh.size):
    area[:,i] = area_t

color_t = rh
colors = np.zeros((vol.size,rh.size),float)
for i in np.arange(vol.size):
    colors[i,:] = color_t

fig = plt.figure(figsize=(6,9))

#fig,axes = plt.subplots(nrows=3,ncols=2)

ax1=plt.subplot(3,2,1)
ax1.scatter(mec_2im, mec_3im, c=colors, cmap='RdYlBu', s=area, alpha=0.8,
        linewidths=0.5, marker='o')

ax1.set_xlim([2, 18])
ax1.set_ylim([2, 18])

lims = [
    np.min([plt.xlim(), plt.ylim()]),
    np.max([plt.ylim(), plt.ylim()]),
    ]

ax1.plot(lims,lims,'k-',alpha=0.75,zorder=0)

ax1.set_ylabel('Homo. Mixing: BC+SUL+OC', fontsize=8)
ax1.set_title('Mass Extinction Coefficient ($m^2/g$)', fontsize=10)
ax1.set_xlabel('Homo. Mixing: BC+SUL only', fontsize=8)
ax1.set_aspect('equal')

plt.grid(True)
plt.tick_params(labelsize=7)


plt.subplot(3,2,2)
plt.scatter(mec_2cs, mec_3cs, c=colors, cmap='RdYlBu', s=area, alpha=0.8, linewidths=0.5, marker='v')
plt.ylabel('Core-shell Mixing: BC+SUL+OC', fontsize=8)
plt.xlabel('Core-shel Mixing: BC+SUL only', fontsize=8)
plt.xlim([2, 18])
plt.ylim([2, 18])
lims = [
    np.min([plt.xlim(), plt.ylim()]),
    np.max([plt.ylim(), plt.ylim()]),
    ]
plt.plot(lims,lims,'k-',alpha=0.75,zorder=0)
plt.grid(True)
plt.tick_params(labelsize=7)
ax2 = plt.subplot(3,2,2)
ax2.set_aspect('equal')

plt.subplot(3,2,3)
plt.scatter(ssa_2im, ssa_3im, c=colors, cmap='RdYlBu', s=area, alpha=0.8, linewidths=0.5, marker='o')
plt.xlim([.3, 1.])
plt.ylim([.3, 1.])
lims = [
    np.min([plt.xlim(), plt.ylim()]),
    np.max([plt.ylim(), plt.ylim()]),
    ]
plt.plot(lims,lims,'k-',alpha=0.75,zorder=0)
plt.ylabel('Homo. Mixing: BC+SUL+OC', fontsize=8)
plt.xlabel('Homo. Mixing: BC+SUL only', fontsize=8)
plt.title('Single Scattering Albedo', fontsize=10)
plt.grid(True)
plt.tick_params(labelsize=7)
ax3 = plt.subplot(3,2,3)
ax3.set_aspect('equal')


plt.subplot(3,2,4)
plt.scatter(ssa_2cs, ssa_3cs, c=colors, cmap='RdYlBu', s=area, alpha=0.8,
        linewidths=0.5, marker='v')
plt.xlim([.3, 1.])
plt.ylim([.3, 1.])
lims = [
    np.min([plt.xlim(), plt.ylim()]),
    np.max([plt.ylim(), plt.ylim()]),
    ]
plt.plot(lims,lims,'k-',alpha=0.75,zorder=0)
plt.ylabel('Core-shell Mixing: BC+SUL+OC', fontsize=8)
plt.xlabel('Core-shell Mixing: BC+SUL only', fontsize=8)
plt.grid(True)
plt.tick_params(labelsize=7)
ax4 = plt.subplot(3,2,4)
ax4.set_aspect('equal')



plt.subplot(3,2,5)
plt.scatter(g_2im, g_3im, c=colors, cmap='RdYlBu', s=area, alpha=0.8,
        linewidths=0.5, marker='o')
plt.xlim([0.35, 0.85])
plt.ylim([0.35, 0.85])
lims = [
    np.min([plt.xlim(), plt.ylim()]),
    np.max([plt.ylim(), plt.ylim()]),
    ]
plt.plot(lims,lims,'k-',alpha=0.75,zorder=0)
plt.ylabel('Homo. Mixing: BC+SUL+OC', fontsize=8)
plt.xlabel('Homo. Mixing: BC+SUL only', fontsize=8)
plt.title('Asymmetry Factor', fontsize=10)
plt.grid(True)
plt.tick_params(labelsize=7)
ax5 = plt.subplot(3,2,5)
ax5.set_aspect('equal')


plt.subplot(3,2,6)
im = plt.scatter(g_2cs, g_3cs, c=colors, cmap='RdYlBu', s=area, alpha=0.8,
        linewidths=0.5, marker='v')
plt.xlim([0.35, 0.85])
plt.ylim([0.35, 0.85])
lims = [
    np.min([plt.xlim(), plt.ylim()]),
    np.max([plt.ylim(), plt.ylim()]),
    ]
plt.tick_params(labelsize=7)
plt.plot(lims,lims,'k-',alpha=0.75,zorder=0)
plt.ylabel('Core-shell Mixing: BC+SUL+OC', fontsize=8)
plt.xlabel('Core-shell Mixing: BC+SUL only', fontsize=8)
plt.grid(True)
ax6 = plt.subplot(3,2,6)
ax6.set_aspect('equal')

plt.tight_layout()

fig.subplots_adjust(bottom=0.13)
cax = fig.add_axes([0.13, 0.05, 0.8, 0.015])
cbar = fig.colorbar(im,cax=cax,ticks=rhpkd,orientation='horizontal')
cbar.solids.set_edgecolor("face")
cbar.set_label("RH (%)")


plt.savefig('radppt_3spec_'+str(volspkd)+'%sul.pdf',format='pdf')
#plt.show()
