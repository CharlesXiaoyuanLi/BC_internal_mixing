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
sup1 = '-ndist'
sup  = '-0.0118'

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
print wlpkdind

# Read in data
mac_im_orig = fh_im.variables['beta_e'][:][:][:] * (1. -
        np.array(fh_im.variables['ssa'][:][:][:])) 
mac_em_orig = fh_em.variables['beta_e'][:][:][:] * (1. -
        np.array(fh_em.variables['ssa'][:][:][:]))
mac_cs_orig = fh_cs.variables['beta_e'][:][:][:] * (1. -
        np.array(fh_cs.variables['ssa'][:][:][:]))

# Calculate the index
rhpkd = np.array([30,60,80,90,94])
rhpkdind = [np.where( rh_orig == s )[0].item() for s in rhpkd]

# Assign picked values
volpkd = np.array([70,80,84,90,94,98])
volpkdind = [np.where( vol_orig == s )[0].item() for s in volpkd]
volcpkdind = [np.where( volc_orig == s )[0].item() for s in volpkd]

print 'vol = '+str(vol_orig[volpkdind])
print 'volc = '+str(volc_orig[volcpkdind])
print 'rh = '+str(rh_orig[rhpkdind])

mac_im_t = mac_im_orig[wlpkdind,volpkdind,:]
mac_im = mac_im_t[:,rhpkdind]
mac_em_t = mac_em_orig[wlpkdind,volpkdind,:]
mac_em = mac_em_t[:,rhpkdind]
mac_cs_t = mac_cs_orig[wlpkdind,volcpkdind,:]
mac_cs = mac_cs_t[:,rhpkdind]
print 'shape = '+str(mac_im.shape)
print mac_cs.shape

units = fh_im.variables['beta_e'].units

print bctyp.upper()+' Pure BC IM mac = '+str(mac_im[0,0])+' '+units
print vol_orig[-1]
print rh_orig[0]
print bctyp.upper()+' Pure BC SSA = '+str(fh_im.variables['ssa'][wlpkdind][-1][0])

fh_im.close()
fh_em.close()
fh_cs.close()

print mac_im[:,0]

# Plot

import matplotlib.pyplot as plt
import matplotlib.lines as lines

vol = vol_orig[volpkdind]
rh = rh_orig[rhpkdind]

area_t = 15 * (vol/25.)**2
area = np.zeros((vol.size,rh.size),float)
for i in np.arange(rh.size):
    area[:,i] = area_t
print area_t

color_t = rh
colors = np.zeros((vol.size,rh.size),float)
for i in np.arange(vol.size):
    colors[i,:] = color_t


#fig, ax = plt.subplots()

fig = plt.figure() #figsize=(7,7)

ax = fig.add_subplot(111)

plt.plot([0, 3.5],[3.5, 3.5], linewidth=.8, alpha=0.2, color='black',
        linestyle=':')

im=plt.scatter(mac_im, mac_cs, c=colors, cmap='RdYlBu', s=area, alpha=0.8,
        linewidths=0.5, marker='o')
#plt.scatter(mac_em, mac_cs, c=colors, cmap='RdYlBu', s=area, alpha=0.8,
#        linewidths=0.5, marker='^')

lims = [
    np.min([plt.xlim(), plt.ylim()]),
    np.max([plt.ylim(), plt.ylim()]),
    ]

plt.plot(lims,lims,'k-',alpha=0.75,zorder=0)

#extraticks = [3.5]
#plt.yticks(list(plt.yticks()[0]) + extraticks)

plt.xlim([0, 3.5])
plt.ylim([0, 3.5])

plt.ylabel('Core-shell Mixing', fontsize=15)
plt.xlabel('Homogeneous Mixing', fontsize=15)
plt.title('BC + Sulfate  Mass Absorption Coefficient ($m^2/g$)')

plt.grid(True)

fig.subplots_adjust(right=0.85)
cax = fig.add_axes([0.88, 0.1, 0.03, 0.8])
cbar = fig.colorbar(im,cax=cax,ticks=rhpkd)
cbar.solids.set_edgecolor("face")
cbar.set_label("RH (%)")


plt.savefig('mac_'+bctyp+'_cs-hm'+sup1+'.pdf',format='pdf')
plt.show()
