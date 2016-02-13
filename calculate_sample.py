from netCDF4 import Dataset
import numpy as np
import scipy.ndimage
import types

### Read data files

# Open netCDF file
nc_file = '../macc_aerosol_rf_201312.nc'
fh = Dataset(nc_file, 'r')

time = fh.variables['t'][:]
lonst = fh.variables['longitude'][:]
latst = fh.variables['latitude'][:]
anthsfct = fh.variables['anthsrf'][:][:][:]

ntime = time.size

nlat = latst.size
nlon = lonst.size

# Extract data from the following domain
lonllim = 70
lonrlim = 140
latblim = 10
latulim = 54

# Calculate the index boundaries of the domain
lonlind = np.where((lonst >= lonllim) & (lonst < lonllim+1))[0]
lonrind = np.where((lonst >lonrlim-1) & (lonst <= lonrlim))[0]
latbind = np.where((latst >= latblim) & (latst < latblim+1))[0]
latuind = np.where((latst > latulim-1) & (latst <= latulim))[0]

lons = lonst[lonlind:lonrind]
lats = latst[latuind:latbind]
anthsfc = anthsfct[:,latuind:latbind,lonlind:lonrind]

print 'Domain: longitude('+str(lonllim)+'~'+str(lonrlim)+'),'+'latitude('+str(latblim)+'~'+str(latulim)+')'
print 'nlon = '+str(lons.size)
print 'nlat = '+str(lats.size)

units = fh.variables['anthsrf'].units

fh.close()

### Calculate total power at each grid

# Calcualte the area weight array
R = 6378100  # radius of the earth in [m] http://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html 
# aw = length in lat  *  length in lon
aw = 2. * np.pi * R * np.cos(np.pi / 180. * lats) / nlon    *    np.pi / 180. * (lats[0] - lats[1]) * R

# Convert W m-2 to kW at each grid point
anthsfc_power = np.zeros( anthsfc.shape )
for i in np.arange(ntime):
    temp = aw * anthsfc[i,:,:].T / 1000.
    anthsfc_power[i,:,:] = temp.T


### Create China Provincial Matrix

# Create a 2D array representing lat x lon
cn_pro_map_t = np.zeros( (nlat,nlon) )

# Read value into the 2D array from the file
f = open("../cn_province_"+str(nlat)+"x"+str(nlon)+".txt","r")
fstring = f.read().split()
istr = 0
for i in range(nlat):
    for j in range(nlon):
        cn_pro_map_t[i,j] = int(fstring[istr])
        istr += 1

# Only take the part for the required domain
cn_pro_map = cn_pro_map_t[latuind:latbind,lonlind:lonrind]

f.close()


### Read in China Province Namelist

f = open("../cn_province_nmlist.txt","r")
pro_nm = f.read().split("\n")
f.close()

print pro_nm

### Calculate total area and power of each province

pro_area = np.zeros(31)
pro_power = np.zeros((ntime,31))

drop = np.zeros(31)

for i in np.arange(cn_pro_map.shape[0]):
    for j in np.arange(cn_pro_map.shape[1]):
        
        if cn_pro_map[i,j] != 0:
            aind = cn_pro_map[i,j]

            pro_area[aind-1] += aw[i]  # total area
            pro_power[:,aind-1] += anthsfc_power[:,i,j]  # total power

            drop[aind-1] += 1

print drop
print pro_area
print pro_power.shape


### 




# Plot

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

m = Basemap(projection='cyl', llcrnrlat=lats[-1], urcrnrlat=lats[0],\
        llcrnrlon=lons[0], urcrnrlon=lons[-1], resolution='c')
# latitudes in the data are in reverse order


lon, lat = np.meshgrid(lons,lats)

# cs = m.contourf(lon,lat,anthsrf,levels=np.arange(-70,5,5))

cs = m.contourf(lon,lat,cn_pro_map,levels=np.arange(0,33,1))

m.drawcoastlines()
m.drawcountries()
m.drawparallels(np.arange(-80,81,20),labels=[1,1,0,0])
m.drawmeridians(np.arange(0,360,30),labels=[0,0,0,1])


cbar = m.colorbar(cs, location='bottom', pad="10%")
cbar.set_label(units)

plt.title('Surface Forcing due to anthropogenic aerosols')
plt.savefig('test_plot',format='ps')
plt.show()
