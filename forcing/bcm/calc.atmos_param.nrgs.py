#!/bin/python2

### Calculate surface albedo

import numpy as np
from netCDF4 import Dataset
import sys


years = [1860,1890,1910,1930,1950,1970,1990]

nyr = len(years)

alb_sfc_nrgs    = np.zeros((15,nyr,12))
f0_nrgs         = np.zeros((15,nyr,12))
ta_nrgs         = np.zeros((15,nyr,12))
rh_nrgs         = np.zeros((15,nyr,12))
### Define constants

re  = 6.37122E+6
d2r = np.pi/180.


for year in years:

    yrn = years.index(year)

    ### Read data

    nc_file = "../albedo/"+str(year)+"_4albedo.nc"

    fh = Dataset(nc_file, 'r')

    f0      = fh.variables['swdn_toa_clr'][:][:][:]

    fsd = fh.variables['swdn_sfc_clr'][:][:][:]
    fsu = fh.variables['swup_sfc_clr'][:][:][:]

    xs,ys,zs = np.where(fsd == 0)
    for x,y,z in zip(xs,ys,zs):
        fsd[x,y,z] = 1

    alb_sfc = np.divide(fsu,fsd)*100.


    lat  = fh.variables['lat'][:]
    latb = fh.variables['latb'][:]
    lon  = fh.variables['lon'][:]

    nlat = len(lat)
    nlon = len(lon)

    time = fh.variables['time'][:]

    fh.close()

    f0_ta      = f0
    fs_ta      = fsd

    ### Read RH data

    nc_file = "../rh/"+str(year)+"_4rh.nc"

    fh = Dataset(nc_file, 'r')

    rh2o = fh.variables['rh2o'][:][:][:][:]   # (time,pfull,lat,lon)
    temp = fh.variables['temp'][:][:][:][:]
    
    pfull = fh.variables['pfull'][:]
    phalf = fh.variables['phalf'][:]

    npres = len(pfull)  # phalf describes the upper and lower boundary
                            # pressure of the grid box.
    fh.close()

    ppkdind = np.where( pfull > 900 )[0]
    pres = pfull[ ppkdind ]
    pres_half = phalf[ np.append(ppkdind, ppkdind[-1]+1) ]
    dpres = np.array([ pres_half[i+1] - pres_half[i] for i in range(len(ppkdind))
        ])
    ### Calculate area of each grid
    weight = np.sin(latb[1:nlat+1]*d2r) - np.sin(latb[0:nlat]*d2r)
    area = 2. * np.pi * np.power(re,2) * weight / nlon
    area_globe = np.array([ area for x in range(nlon) ]).transpose()

    regions = ["Global","Eastern Asia","Southeastern Asia","Southern Asia", \
            "Northern America","Central America","South America","Europe", \
            "Northern Africa","Western Africa","Eastern Africa","Southern Africa",\
            "Middle Africa","Pacific Warm Pool","Arctic"]

    regionids = [0, 8, 19, 21, 15, 5, 18, 25, 14, 23, 7, 20, 13, 26, 27] 

    ### Read mask ncfile
    mask_file = "world_regions_2008_"+str(nlat)+"x"+str(nlon)+".nc"
    fh_mask = Dataset(mask_file, 'r')
    rgmask = fh_mask.variables['mask'][:][:]

    fh_mask.close()

    ### Calculate average for each region
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

        rg_alb_sfc_array = np.array([ alb_sfc[:,x,y] for x,y in zip(ilats,ilons) ])
        rg_f0_array = np.array([ f0[:,x,y] for x,y in zip(ilats,ilons) ])

        rg_f0_ta_array = np.array([ f0_ta[:,x,y] for x,y in zip(ilats,ilons) ])
        rg_fs_ta_array = np.array([ fs_ta[:,x,y] for x,y in zip(ilats,ilons) ])

        rg_rh2o_array = np.array([ rh2o[:,:,x,y] for x,y in zip(ilats,ilons) ])
        rg_temp_array = np.array([ temp[:,:,x,y] for x,y in zip(ilats,ilons) ])

        rg_area_array = np.array([ area_globe[x,y] for x,y in zip(ilats,ilons) ])

        rg_alb_sfc = np.average(rg_alb_sfc_array, axis=0, weights=rg_area_array)
        rg_f0 = np.average(rg_f0_array, axis=0, weights=rg_area_array)
        rg_f0_ta = np.average(rg_f0_ta_array, axis=0, weights=rg_area_array)
        rg_fs_ta = np.average(rg_fs_ta_array, axis=0, weights=rg_area_array) # should be from nobcsul!
        rg_ta = rg_fs_ta / rg_f0_ta
        
        rg_rh2o_pres = np.average(rg_rh2o_array, axis=0, weights=rg_area_array)[:,ppkdind]
        rg_temp_pres = np.average(rg_temp_array, axis=0, weights=rg_area_array)[:,ppkdind]
        
        # Calculate RH for each layer
        # convert mixing ratio to specific humidity
        rg_qv_pres = rg_rh2o_pres / (rg_rh2o_pres + 1)
        # calculate vapor pressure
        rg_e_pres = np.array([ rg_qv_pres[i,:] * pres for i in
            range(rg_qv_pres.shape[0]) ])  # hPa
        # calculate saturated vapor pressure using Clausius-Clapeyron Equation
        rg_es_pres = 6.11 * np.exp(5420 * (1./273 - 1./rg_temp_pres))
        rg_rh_pres = rg_e_pres / rg_es_pres * 100

        # Calculate average RH for the entire bottom atmosphere (<800hPa)
        ## calculate average mixing ratio and temperature
        #rg_rh2o = np.average(rg_rh2o_pres, axis=1, weights=dpres)
        #rg_temp = np.average(rg_temp_pres, axis=1, weights=dpres)
        ## convert mixing ratio to specific humidity
        #rg_qv = rg_rh2o / (rg_rh2o + 1)
        ## calculate vapor pressure
        #rg_e = rg_qv * np.average(pres_half[[0,-1]])
        ## calculate saturated vapor pressure using C-C Equation
        
        # Calculate RH for the bottom layer
        rg_rh2o = rg_rh2o_pres[:,-1]
        rg_temp = rg_temp_pres[:,-1]
        rg_qv = rg_rh2o / (rg_rh2o + 1)
        rg_e = rg_qv * pres[-1]
       
        rg_es = 6.11 * np.exp(5420 * (1./273 - 1./rg_temp))
        rg_rh = rg_e / rg_es * 100

        print region
        rgind = regions.index(region)
        alb_sfc_nrgs[rgind,yrn,:] = rg_alb_sfc
        f0_nrgs[rgind,yrn,:] = rg_f0
        ta_nrgs[rgind,yrn,:] = rg_ta
        rh_nrgs[rgind,yrn,:] = rg_rh

### Write data to ncfile
outfile = Dataset('atmos_param_'+str(len(regions))+'rgs.nc','w',format='NETCDF4')

outfile.createDimension('region',len(regions))
outfile.createDimension('year',nyr)
months = np.arange(12) + 1
outfile.createDimension('month',len(months))

year_dim  = outfile.createVariable('year',np.int32,('year',))
month_dim = outfile.createVariable('month',np.int32,('month',))

region_dim = outfile.createVariable('region',np.int32,('region',))
description = ''
for region, rgid in zip(regions, regionids):
    description += region + ': ' + str(rgid) + ',  '

region_dim.description =description

alb_sfc         = outfile.createVariable('alb_sfc',np.float32,('region','year','month'))
alb_sfc.long_name = 'surface albedo'
alb_sfc.units   = '%'

f0              = outfile.createVariable('f0',np.float32,('region','year','month'))
f0.long_name    = 'SW flux down at TOA'
f0.units       = 'watts/m2'

ta              = outfile.createVariable('ta',np.float32,('region','year','month'))
ta.long_name    = 'atmospheric transmittance'

rh              = outfile.createVariable('rh',np.float32,('region','year','month'))
rh.long_name    = 'relative humidity'
rh.units        = '%'

region_dim[:]   = regionids
year_dim[:]     = years
month_dim[:]    = months

alb_sfc[:]  = alb_sfc_nrgs
f0[:]       = f0_nrgs
ta[:]       = ta_nrgs
rh[:]       = rh_nrgs

outfile.close()


