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
#ncfile_im = 'sul-'+bctyp+'-im.nc'
ncfile_im = 'sul_bc_im.nc'
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
mec_im_orig = fh_im.variables['beta_e'][:][:][:]
mec_em_orig = fh_em.variables['beta_e'][:][:][:]
mec_cs_orig = fh_cs.variables['beta_e'][:][:][:]

# Calculate the index
rhpkd = np.array([30,40,50,60,70,80,84,88,91,93,96])
rhpkdind = [np.where( rh_orig == s )[0].item() for s in rhpkd]

# Assign picked values
wl  = wlpkd
volint = -3
vol = vol_orig[::volint]
print 'vol = '+str(vol)
rh  = rhpkd
print 'rh = '+str(rh)

mec_im_t = mec_im_orig[wlpkdind,::volint,:]
mec_im = mec_im_t[:,rhpkdind]
mec_em_t = mec_em_orig[wlpkdind,::volint,:]
mec_em = mec_em_t[:,rhpkdind]
mec_cs_t = mec_cs_orig[wlpkdind,::volint,:]
mec_cs = mec_cs_t[:,rhpkdind]
print 'shape = '+str(mec_im.shape)

units = fh_im.variables['beta_e'].units

print bctyp.upper()+' Pure BC MEC = '+str(mec_im[0,0])+' '+units

fh_im.close()
fh_em.close()
fh_cs.close()

# Plot

import matplotlib.pyplot as plt

area_t = 15 * ((100-vol)/25)**2
area = np.zeros((vol.size,rh.size),float)
for i in np.arange(rh.size):
    area[:,i] = area_t

color_t = rh
colors = np.zeros((vol.size,rh.size),float)
for i in np.arange(vol.size):
    colors[i,:] = color_t


fig, ax = plt.subplots()

plt.scatter(mec_em, mec_im, c=colors, cmap='RdYlBu', s=area, alpha=0.8,
        linewidths=0.5, marker='o')
#plt.scatter(mec_em, mec_cs, c=colors, cmap='RdYlBu', s=area, alpha=0.8,
#        linewidths=0.5, marker='^')

plt.xlim([3, 14])
plt.ylim([3, 14])

lims = [
    np.min([plt.xlim(), plt.ylim()]),
    np.max([plt.ylim(), plt.ylim()]),
    ]

plt.plot(lims,lims,'k-',alpha=0.75,zorder=0)

plt.xlabel('External Mixing', fontsize=15)
plt.ylabel('Homogeneous Internal Mixing', fontsize=15)
plt.title('BC + Sulfate  Mass Extinction Coefficient ($m^2/g$)')

plt.grid(True)

plt.savefig('mec_'+bctyp+'.pdf',format='pdf')
plt.show()
