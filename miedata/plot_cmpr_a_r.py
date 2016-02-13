#!/bin/python2

### This ncfile is set up to compare the absorptivity and reflectivity of
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

gr = 1.

sup = '-0.0118'
sup1 = '-ndist'

# Open netCDF file
ncfile_im = 'sul-'+bctyp+'-im'+sup1+'.nc'
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
print wlpkdind

# Read in data
a_im_orig = fh_im.variables['beta_e'][:][:][:] * (1. -
        np.array(fh_im.variables['ssa'][:][:][:])) 
a_em_orig = fh_em.variables['beta_e'][:][:][:] * (1. -
        np.array(fh_em.variables['ssa'][:][:][:]))
a_cs_orig = fh_cs.variables['beta_e'][:][:][:] * (1. -
        np.array(fh_cs.variables['ssa'][:][:][:]))

r_im_orig = fh_im.variables['beta_e'][:][:][:] * fh_im.variables['ssa'][:][:][:]\
* 0.5 * (1. - gr * fh_im.variables['g'][:][:][:])
r_em_orig = fh_em.variables['beta_e'][:][:][:] * fh_em.variables['ssa'][:][:][:]\
* 0.5 * (1. - gr * fh_em.variables['g'][:][:][:])
r_cs_orig = fh_cs.variables['beta_e'][:][:][:] * fh_cs.variables['ssa'][:][:][:]\
* 0.5 * (1. - gr * fh_cs.variables['g'][:][:][:])


# Calculate the index
rhpkd = np.array([30,60,80,90,94])
rhpkdind = [np.where( rh_orig == s )[0].item() for s in rhpkd]

# Assign picked values
volpkd = np.array([70,75,80,84,88,92,96,98])
volpkdind = [np.where( vol_orig == s )[0].item() for s in volpkd]
volcpkdind = [np.where( volc_orig == s )[0].item() for s in volpkd]

print 'vol = '+str(vol_orig[volpkdind])
print 'volc = '+str(volc_orig[volcpkdind])
print 'rh = '+str(rh_orig[rhpkdind])

a_im_t = a_im_orig[wlpkdind,volpkdind,:]
a_im = a_im_t[:,rhpkdind]
a_em_t = a_em_orig[wlpkdind,volpkdind,:]
a_em = a_em_t[:,rhpkdind]
a_cs_t = a_cs_orig[wlpkdind,volcpkdind,:]
a_cs = a_cs_t[:,rhpkdind]

r_im_t = r_im_orig[wlpkdind,volpkdind,:]
r_im = r_im_t[:,rhpkdind]
r_em_t = r_em_orig[wlpkdind,volpkdind,:]
r_em = r_em_t[:,rhpkdind]
r_cs_t = r_cs_orig[wlpkdind,volcpkdind,:]
r_cs = r_cs_t[:,rhpkdind]

da_im = a_im - a_em
dr_im = r_em - r_im

da_cs = a_cs - a_em
dr_cs = r_em - r_cs

print 'shape = '+str(a_im.shape)
print a_cs.shape

units = fh_im.variables['beta_e'].units

fh_im.close()
fh_em.close()
fh_cs.close()

print a_im[:,0]

# Plot

import matplotlib.pyplot as plt

vol = vol_orig[volpkdind]
rh = rh_orig[rhpkdind]

area_t = 15 * (vol/30.)**2
area = np.zeros((vol.size,rh.size),float)
for i in np.arange(rh.size):
    area[:,i] = area_t
print area_t

color_t = rh
colors = np.zeros((vol.size,rh.size),float)
for i in np.arange(vol.size):
    colors[i,:] = color_t

fig, axes = plt.subplots(nrows=1,ncols=2)

plt.subplot(1,2,1)

plt.scatter(dr_im, da_im, c=colors, cmap='RdYlBu', s=area, alpha=0.8,
        linewidths=0.5, marker='o')

#plt.xlim([-0.4, 0.3])
#plt.ylim([0, 1.6]) 

#plt.plot([0, 0.8],[0, 3.5],'k-',alpha=0.5,zorder=0)

plt.ylabel('$\Delta a$', fontsize=20)
plt.xlabel('$- \Delta r$', fontsize=20)
plt.title('Homogeneous Mixing',loc='left')
plt.grid(True)

plt.subplot(1,2,2)

im = plt.scatter(dr_cs, da_cs, c=colors, cmap='RdYlBu', s=area, alpha=0.8,
        linewidths=0.5, marker='^')

#plt.xlim([-0.4, 0.3])
#plt.ylim([0, 1.6]) 

lims = [
    np.min([plt.xlim(), plt.ylim()]),
    np.max([plt.ylim(), plt.ylim()]),
    ]
#plt.plot(lims,lims,'k-',alpha=0.75,zorder=0)

plt.ylabel('$\Delta a$', fontsize=20)
plt.xlabel('$- \Delta r$', fontsize=20)
plt.title('Core-shell Mixing',loc='left')
plt.grid(True)

plt.tight_layout()

fig.subplots_adjust(right=0.85)
cax = fig.add_axes([0.9, 0.13, 0.02, 0.8])
cbar = fig.colorbar(im,cax=cax,ticks=rhpkd,orientation='vertical')
cbar.solids.set_edgecolor("face")
cbar.set_label("RH (%)")


#plt.savefig('dadr_'+bctyp+'_im'+sup1+'.pdf',format='pdf')
#For magnified figure
plt.savefig('dadr_'+bctyp+'_im'+sup1+'.pdf',format='pdf')

plt.show()
