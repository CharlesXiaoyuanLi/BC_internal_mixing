from netCDF4 import Dataset
import numpy as np
import scipy.ndimage
import types

# Open netCDF file
nc_file = 'macc_aerosol_rf_201312.nc'
fh = Dataset(nc_file, 'r')

lonst = fh.variables['longitude'][:]
latst = fh.variables['latitude'][:]
anthsfct = fh.variables['anthsrf'][23][:][:]

lonllim = 70
lonrlim = 140
latblim = 10
latulim = 54


lonlind = np.where((lonst >= lonllim) & (lonst < lonllim+1))[0]
lonrind = np.where((lonst >lonrlim-1) & (lonst <= lonrlim))[0]
latbind = np.where((latst >= latblim) & (latst < latblim+1))[0]
latuind = np.where((latst > latulim-1) & (latst <= latulim))[0]

lons = lonst[lonlind:lonrind]
lats = latst[latuind:latbind]
anthsfc = anthsfct[latuind:latbind,lonlind:lonrind]

print 'Domain: longitude('+str(lonllim)+'~'+str(lonrlim)+'),'+'latitude('+str(latblim)+'~'+str(latulim)+')'
print 'nlon = '+str(len(lons))
print 'nlat = '+str(len(lats))

units = fh.variables['anthsrf'].units

fh.close()

# Read the China Provincial Matrix

# Create a 2D array representing lat x lon
nlat = 161
nlon = 320
#cn_pro_map = [[0 for x in range(nlon)] for x in range(nlat)]
cn_pro_map = anthsfct
f = open("cn_province_"+str(nlat)+"x"+str(nlon)+".txt","r")

fstring = f.read().split()
istr = 0
for i in range(nlat):
    for j in range(nlon):
        cn_pro_map[i,j] = int(fstring[istr])
        istr += 1

print cn_pro_map

# Plot

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

m = Basemap(projection='cyl', llcrnrlat=lats[-1], urcrnrlat=lats[0],\
        llcrnrlon=lons[0], urcrnrlon=lons[-1], resolution='c')
# latitudes in the data are in reverse order


lon, lat = np.meshgrid(lons,lats)

# cs = m.contourf(lon,lat,anthsrf,levels=np.arange(-70,5,5))

cs = m.contourf(lon,lat,cn_pro_map[latuind:latbind,lonlind:lonrind],levels=np.arange(0,33,1))

m.drawcoastlines()
m.drawcountries()
m.drawparallels(np.arange(-80,81,20),labels=[1,1,0,0])
m.drawmeridians(np.arange(0,360,30),labels=[0,0,0,1])


cbar = m.colorbar(cs, location='bottom', pad="10%")
cbar.set_label(units)

plt.title('Surface Forcing due to anthropogenic aerosols')
plt.savefig('test_plot',format='ps')
plt.show()
