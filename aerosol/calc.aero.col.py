### Calculate aerosol column density

import numpy as np
from netCDF4 import Dataset

### Define constants
re  = 6.37122E+6
d2r = np.pi/180.


### Read data from ncfile
fn = 'aerosol.used.nc'
fh = Dataset(fn, 'r')

bc_raw  = fh.variables['black_carbon'][:][:][:][:]
so4_raw = fh.variables['so4_anthro'][:][:][:][:]
oc_raw  = fh.variables['organic_carbon'][:][:][:][:]

time = fh.variables['time'][:]
lat  = fh.variables['lat'][:]
lon  = fh.variables['lon'][:]
latb = fh.variables['latb'][:]

ntime = len(time)
nlat  = len(lat)
nlon  = len(lon)

fh.close()

# Add up all vertical levels together for column density
bc_col_raw  = np.sum(bc_raw, axis=1)
so4_col_raw = np.sum(so4_raw, axis=1)
oc_col_raw  = np.sum(oc_raw, axis=1)
print bc_col_raw.shape

# separate year and month into two coordinates
bc_col  = np.array([ bc_col_raw[i:i+12,:,:] for i in np.arange(0,ntime,12) ])
so4_col = np.array([ so4_col_raw[i:i+12,:,:] for i in np.arange(0,ntime,12) ])
oc_col  = np.array([ oc_col_raw[i:i+12,:,:] for i in np.arange(0,ntime,12) ])
print bc_col.shape

### Calculate mean column density for each region
# Calculate area of each grid
weight = np.sin(latb[1:nlat+1]*d2r) - np.sin(latb[0:nlat]*d2r)
area = 2. * np.pi * np.power(re,2) * weight / nlon
area_globe = np.array([ area for x in range(nlon) ]).transpose()

regions = ["Global","Eastern Asia","Southeastern Asia","Southern Asia", \
        "Northern America","Central America","South America","Europe", \
        "Northern Africa","Western Africa","Eastern Africa","Southern Africa",\
        "Middle Africa","Pacific Warm Pool","Arctic"]

regionids = [0, 8, 19, 21, 15, 5, 18, 25, 14, 23, 7, 20, 13, 26, 27] 

# Read mask ncfile
mask_file = "../forcing/bcm/world_regions_2008_"+str(nlat)+"x"+str(nlon)+".nc"
fh_mask = Dataset(mask_file, 'r')
rgmask = fh_mask.variables['mask'][:][:]

fh_mask.close()

# Calculate average for each region
bc_col_nrgs  = np.zeros((15,bc_col.shape[0],bc_col.shape[1]))
so4_col_nrgs = np.zeros((15,so4_col.shape[0],so4_col.shape[1]))
oc_col_nrgs  = np.zeros((15,oc_col.shape[0],oc_col.shape[1]))

for region, rgid in zip(regions, regionids):

    if ( region == 'Global' ):
        ilats = [ x for x in range(nlat) for y in range(nlon) ]
        ilons = [ y for x in range(nlat) for y in range(nlon) ]
    elif ( region == 'Pacific Warm Pool' ):
        lats = np.where( (lat >= -10) & (lat <= 10) )[0]
        lons = np.where( (lon >= 90) & (lon <= 130) )[0]
        ilats = [ x for x in lats for y in lons ]
        ilons = [ y for x in lats for y in lons ]
    elif ( region == 'Arctic' ):
        lats = np.where( lat > 66 )[0]
        lons = range(nlon)
        ilats = [ x for x in lats for y in lons ]
        ilons = [ y for x in lats for y in lons ]
    elif ( region == 'Europe' ):
        ilats, ilons = np.where( (rgmask == 25) | (rgmask == 22) | (rgmask == \
            16) | (rgmask == 9) )
    else:
        ilats, ilons = np.where( rgmask == rgid )

    rg_bc_col_array = np.array([ bc_col[:,:,x,y] for x,y in zip(ilats,ilons) ])
    rg_so4_col_array = np.array([ so4_col[:,:,x,y] for x,y in zip(ilats,ilons) ])
    rg_oc_col_array = np.array([ oc_col[:,:,x,y] for x,y in zip(ilats,ilons) ])

    rg_area_array = np.array([ area_globe[x,y] for x,y in zip(ilats,ilons) ])

    rg_bc_col = np.average(rg_bc_col_array, axis=0, weights=rg_area_array)
    rg_so4_col = np.average(rg_so4_col_array, axis=0, weights=rg_area_array)
    rg_oc_col = np.average(rg_oc_col_array, axis=0, weights=rg_area_array)

    rgind = regions.index(region)
    bc_col_nrgs[rgind,:,:] = rg_bc_col
    so4_col_nrgs[rgind,:,:] = rg_so4_col
    oc_col_nrgs[rgind,:,:] = rg_oc_col

# Print column density in 1990 as an example

for region in regions:

    rgind = regions.index(region)
    
    print region+': '
    print 'BC Column Density: ',
    print bc_col_nrgs[rgind,-2,:]*1E9
    print 'Sulfate Column Density: ',
    print so4_col_nrgs[rgind,-2,:]*1E9
    print 'Volume ratio of sulfate: ',
    print so4_col_nrgs[rgind,-2,:]/1.74/(bc_col_nrgs[rgind,-2,:]/1.8 +
            so4_col_nrgs[rgind,-2,:]/1.74)
    print '\n'

### Write output to ncfile
outfile = Dataset('aerocol_'+str(len(regions))+'rgs.nc','w',format='NETCDF4')

outfile.createDimension('region',len(regions))
outfile.createDimension('year',8)
outfile.createDimension('month',12)

years  = np.array([1860,1890,1910,1930,1950,1970,1990,2002])
month = range(1,13)

year_dim = outfile.createVariable('year',np.int32,('year',))
month_dim = outfile.createVariable('month',np.int32,('month',))

region_dim = outfile.createVariable('region',np.int32,('region',))
description = ''
for region, rgid in zip(regions, regionids):
    description += region + ': ' + str(rgid) + ',  '

region_dim.description = description

bc_col_nc  = outfile.createVariable('bc_col',np.float32,('region','year','month'))
bc_col_nc.units  = "kg/m2"
so4_col_nc = outfile.createVariable('so4_col',np.float32,('region','year','month'))
so4_col_nc.units = "kg/m2"
oc_col_nc  = outfile.createVariable('oc_col',np.float32,('region','year','month'))
oc_col_nc.units  = "kg/m2"

region_dim[:] = regionids
year_dim[:]   = years
month_dim[:]  = month

bc_col_nc[:]  = bc_col_nrgs
so4_col_nc[:] = so4_col_nrgs
oc_col_nc[:]  = oc_col_nrgs

outfile.close()
