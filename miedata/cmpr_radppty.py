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

sup = '-0.0118'
sup1 = '-ndist'

# Open netCDF file
ncfile_im = 'sul-'+bctyp+'-im'+sup1+'.nc'
#ncfile_im = 'sul_bc_im.nc'
ncfile_em = 'sul-'+bctyp+'-em'+sup+'.nc'
ncfile_cs = 'sul-'+bctyp+'-coat'+sup1+'.new.nc'

fh_im = Dataset(ncfile_im, 'r')
fh_em = Dataset(ncfile_em, 'r')
fh_cs = Dataset(ncfile_cs, 'r')

# Read in coordinates
wl_orig  = fh_im.variables['wl'][:]
vol_orig = fh_im.variables['vol'][:]
rh_orig  = fh_im.variables['RH'][:]

volc_orig = fh_cs.variables['vol'][:]

wlpkd = 0.55  # choose the wavelength at 550nm to represent shortwave radiation
wlpkdind = np.where((wl_orig >= wlpkd-0.02) & (wl_orig <= wlpkd+0.02))[0].item()

# Read in data
mec_im_orig = fh_im.variables['beta_e'][:][:][:]
mec_em_orig = fh_em.variables['beta_e'][:][:][:]
mec_cs_orig = fh_cs.variables['beta_e'][:][:][:]

ssa_im_orig = fh_im.variables['ssa'][:][:][:]
ssa_em_orig = fh_em.variables['ssa'][:][:][:]
ssa_cs_orig = fh_cs.variables['ssa'][:][:][:]

g_im_orig = fh_im.variables['g'][:][:][:]
g_em_orig = fh_em.variables['g'][:][:][:]
g_cs_orig = fh_cs.variables['g'][:][:][:]


# Calculate the index
rhpkd = np.array([30,50,60,70,80,90,94])
rhpkdind = [np.where( rh_orig == s )[0].item() for s in rhpkd]

# Assign picked values
volpkd = np.array([70,80,84,90,94,98])
volpkdind = [np.where( vol_orig == s )[0].item() for s in volpkd]
volcpkdind = [np.where( volc_orig == s )[0].item() for s in volpkd]

print 'vol = '+str(vol_orig[volpkdind])
print 'volc = '+str(volc_orig[volcpkdind])
print 'rh = '+str(rh_orig[rhpkdind])

mec_im_t = mec_im_orig[wlpkdind,volpkdind,:]
mec_im = mec_im_t[:,rhpkdind]
mec_em_t = mec_em_orig[wlpkdind,volpkdind,:]
mec_em = mec_em_t[:,rhpkdind]
mec_cs_t = mec_cs_orig[wlpkdind,volcpkdind,:]
mec_cs = mec_cs_t[:,rhpkdind]

units = fh_im.variables['beta_e'].units

print vol_orig[-1]
print bctyp.upper()+' Pure BC MEC = '+str(mec_em_orig[wlpkdind,-1,0])+' '+units

ssa_im_t = ssa_im_orig[wlpkdind,volpkdind,:]
ssa_im = ssa_im_t[:,rhpkdind]
ssa_em_t = ssa_em_orig[wlpkdind,volpkdind,:]
ssa_em = ssa_em_t[:,rhpkdind]
ssa_cs_t = ssa_cs_orig[wlpkdind,volcpkdind,:]
ssa_cs = ssa_cs_t[:,rhpkdind]

g_im_t = g_im_orig[wlpkdind,volpkdind,:]
g_im = g_im_t[:,rhpkdind]
g_em_t = g_em_orig[wlpkdind,volpkdind,:]
g_em = g_em_t[:,rhpkdind]
g_cs_t = g_cs_orig[wlpkdind,volcpkdind,:]
g_cs = g_cs_t[:,rhpkdind]


fh_im.close()
fh_em.close()
fh_cs.close()

# Plot

import matplotlib.pyplot as plt

vol = vol_orig[volpkdind]
rh  = rh_orig[rhpkdind]

area_t = 15 * (vol/45.)**2
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
ax1.scatter(mec_em, mec_im, c=colors, cmap='RdYlBu', s=area, alpha=0.8,
        linewidths=0.5, marker='o')

#ax1.set_xlim([3, 14])
#ax1.set_ylim([3, 14])
ax1.set_xlim([0, 25])
ax1.set_ylim([0, 25])

lims = [
    np.min([plt.xlim(), plt.ylim()]),
    np.max([plt.ylim(), plt.ylim()]),
    ]

ax1.plot(lims,lims,'k-',alpha=0.75,zorder=0)

ax1.set_ylabel('Homogeneous Mixing', fontsize=8)
ax1.set_title('Mass Extinction Coefficient ($m^2/g$)', fontsize=10)
ax1.set_xlabel('External Mixing', fontsize=8)
ax1.set_aspect('equal')

plt.grid(True)
plt.tick_params(labelsize=7)


plt.subplot(3,2,2)
plt.scatter(mec_em, mec_cs, c=colors, cmap='RdYlBu', s=area, alpha=0.8,
        linewidths=0.5, marker='v')
plt.ylabel('Core-shell Mixing', fontsize=8)
plt.xlabel('External Mixing', fontsize=8)
#plt.xlim([3, 14])
#plt.ylim([3, 14])
plt.xlim([0, 25])
plt.ylim([0, 25])
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
plt.scatter(ssa_em, ssa_im, c=colors, cmap='RdYlBu', s=area, alpha=0.8,
        linewidths=0.5, marker='o')
#plt.xlim([0.35, 1])
#plt.ylim([0.35, 1])
plt.xlim([.3, 1.])
plt.ylim([.3, 1.])
lims = [
    np.min([plt.xlim(), plt.ylim()]),
    np.max([plt.ylim(), plt.ylim()]),
    ]
plt.plot(lims,lims,'k-',alpha=0.75,zorder=0)
plt.ylabel('Homogeneous Mixing', fontsize=8)
plt.xlabel('External Mixing', fontsize=8)
plt.title('Single Scattering Albedo', fontsize=10)
plt.grid(True)
plt.tick_params(labelsize=7)
ax3 = plt.subplot(3,2,3)
ax3.set_aspect('equal')


plt.subplot(3,2,4)
plt.scatter(ssa_em, ssa_cs, c=colors, cmap='RdYlBu', s=area, alpha=0.8,
        linewidths=0.5, marker='v')
#plt.xlim([0.35, 1])
#plt.ylim([0.35, 1])
plt.xlim([.3, 1.])
plt.ylim([.3, 1.])
lims = [
    np.min([plt.xlim(), plt.ylim()]),
    np.max([plt.ylim(), plt.ylim()]),
    ]
plt.plot(lims,lims,'k-',alpha=0.75,zorder=0)
plt.ylabel('Core-shell Mixing', fontsize=8)
plt.xlabel('External Mixing', fontsize=8)
plt.grid(True)
plt.tick_params(labelsize=7)
ax4 = plt.subplot(3,2,4)
ax4.set_aspect('equal')



plt.subplot(3,2,5)
plt.scatter(g_em, g_im, c=colors, cmap='RdYlBu', s=area, alpha=0.8,
        linewidths=0.5, marker='o')
#plt.xlim([0.55, 0.8])
#plt.ylim([0.55, 0.8])
plt.xlim([0.35, 0.85])
plt.ylim([0.35, 0.85])
lims = [
    np.min([plt.xlim(), plt.ylim()]),
    np.max([plt.ylim(), plt.ylim()]),
    ]
plt.plot(lims,lims,'k-',alpha=0.75,zorder=0)
plt.ylabel('Homogeneous Mixing', fontsize=8)
plt.xlabel('External Mixing', fontsize=8)
plt.title('Asymmetry Factor', fontsize=10)
plt.grid(True)
plt.tick_params(labelsize=7)
ax5 = plt.subplot(3,2,5)
ax5.set_aspect('equal')


plt.subplot(3,2,6)
im = plt.scatter(g_em, g_cs, c=colors, cmap='RdYlBu', s=area, alpha=0.8,
        linewidths=0.5, marker='v')
#plt.xlim([0.55, 0.8])
#plt.ylim([0.55, 0.8])
plt.xlim([0.35, 0.85])
plt.ylim([0.35, 0.85])
lims = [
    np.min([plt.xlim(), plt.ylim()]),
    np.max([plt.ylim(), plt.ylim()]),
    ]
plt.tick_params(labelsize=7)
plt.plot(lims,lims,'k-',alpha=0.75,zorder=0)
plt.ylabel('Core-shell Mixing', fontsize=8)
plt.xlabel('External Mixing', fontsize=8)
plt.grid(True)
ax6 = plt.subplot(3,2,6)
ax6.set_aspect('equal')

plt.tight_layout()

fig.subplots_adjust(bottom=0.13)
cax = fig.add_axes([0.13, 0.05, 0.8, 0.015])
cbar = fig.colorbar(im,cax=cax,ticks=rhpkd,orientation='horizontal')
cbar.solids.set_edgecolor("face")
cbar.set_label("RH (%)")


plt.savefig('radppty_'+bctyp+sup+sup1+'.pdf',format='pdf')
plt.show()
