### Calculate the linear regression of TOA vs. SFC

### "Newtheory" refers to different shortwave spectrum between surface and TOA.
#It seems to be not very effective.

import numpy as np
from netCDF4 import Dataset
from scipy import stats
import itertools
import matplotlib.pyplot as plt
import sys

### Read required document names
emtype = raw_input("EMTYPE: Choose between 'embcm' and 'embcm_0.05'\n")
imtype = raw_input("IMTYPE: Choose among 'imbcm', 'imbcm_ndist', 'csbcm' \
and 'csbcm_ndist'\n")
rgs = raw_input("# of regions: 13 or 15?\n")


### Read data of Standalone Model

fn_toa = "toa_"+emtype+"_"+imtype+"_"+rgs+"rgs"+".nc"
fn_sfc = "sfc_"+emtype+"_"+imtype+"_"+rgs+"rgs"+".nc"
fn_par = "atmos_param_"+rgs+"rgs"+".nc"

fh_toa = Dataset(fn_toa, 'r')
fh_sfc = Dataset(fn_sfc, 'r')
fh_par = Dataset(fn_par, 'r')

toa_di = fh_toa.variables['rf_di'][:][:]
toa_int = fh_toa.variables['rf_int'][:][:]
toa_ext = fh_toa.variables['rf_ext'][:][:]
sfc_di = fh_sfc.variables['rf_di'][:][:]
sfc_int = fh_sfc.variables['rf_int'][:][:]
sfc_ext = fh_sfc.variables['rf_ext'][:][:]
f0      = fh_par.variables['f0'][:][:][:]
ta      = fh_par.variables['ta'][:][:][:]
alb_sfc = fh_par.variables['alb_sfc'][:][:][:]
rh_rg   = fh_par.variables['rh'][:][:][:]

fh_toa.close()
fh_sfc.close()
fh_par.close()


### Define regions
regions = ["Global","Eastern Asia","Southeastern Asia","Southern Asia", \
        "Northern America","Central America","South America","Europe", \
        "Northern Africa","Western Africa","Eastern Africa","Southern Africa",\
        "Middle Africa","Pacific Warm Pool","Arctic"]

### Read in column density data
fn_col = "./aerocol_"+rgs+"rgs.nc"
fh_col = Dataset(fn_col, 'r')

bc_col_t  = fh_col.variables['bc_col'][:][:][:] #no 2002 
so4_col_t = fh_col.variables['so4_col'][:][:][:] 

years = fh_col.variables['year'][:]
rgids = fh_col.variables['region'][:]

bc_col = bc_col_t[:,:-1,:]
so4_col = so4_col_t[:,:-1,:]

volr = so4_col/1.74/(so4_col/1.74 + bc_col/1.8)

### Read in MAC & MSC data
ncfile_em = '../../miedata/sul-bcm-em-0.0118.nc'
if (imtype == 'imbcm_ndist'):
    ncfile_cs = '../../miedata/sul-bcm-im-ndist.nc'
elif (imtype == 'csbcm_ndist'):
    ncfile_cs = '../../miedata/sul-bcm-coat-ndist.nc'

fh_em = Dataset(ncfile_em, 'r')
fh_cs = Dataset(ncfile_cs, 'r')

wl      = fh_em.variables['wl'][:]
volem   = fh_em.variables['vol'][:]
rh      = fh_em.variables['RH'][:]
volcs   = fh_cs.variables['vol'][:]

mac_em_t0 = fh_em.variables['beta_e'][:][:][:] * (1. - \
        np.array(fh_em.variables['ssa'][:][:][:]))
mac_cs_t0 = fh_cs.variables['beta_e'][:][:][:] * (1. - \
        np.array(fh_cs.variables['ssa'][:][:][:]))

msc_em_t0 = fh_em.variables['beta_e'][:][:][:] * fh_em.variables['ssa'][:][:][:]
g_em   = fh_em.variables['g'][:][:][:]
msc_cs_t0 = fh_cs.variables['beta_e'][:][:][:] * fh_cs.variables['ssa'][:][:][:]
g_cs   = fh_cs.variables['g'][:][:][:]

mscb_em_t0 = msc_em_t0 * 0.5 * (1. - 0.5*g_em)  # currently not appplied: effective backscattering ratio
mscb_cs_t0 = msc_cs_t0 * 0.5 * (1. - 0.5*g_cs)  # by considering isotropic insolation

mscb_em1_t0 = msc_em_t0 * 0.5 * (1. - 0.5*g_em)  # currently not appplied: effective backscattering ratio
mscb_cs1_t0 = msc_cs_t0 * 0.5 * (1. - 0.5*g_cs)  # by considering isotropic insolation


# Pick the MAC data that match
wlpkdind = np.where( wl < 2.5 )[0]
wlpkdwgt = np.array([0.1*250,0.2*375,0.3*500,0.8*430,1.6*255,1.9*180,1.4*180,0.7*80,0])
wlpkdwgt1 = np.array([0*250,0*375,0.25*500,0.45*430,0.8*255,1.25*180,0.3*180,0*80,0])

print wlpkdind
print wl[wlpkdind]

mac_em_t = np.average(mac_em_t0[wlpkdind,:,:],axis=0,weights=wlpkdwgt)
mac_cs_t = np.average(mac_cs_t0[wlpkdind,:,:],axis=0,weights=wlpkdwgt)
msc_em_t = np.average(msc_em_t0[wlpkdind,:,:],axis=0,weights=wlpkdwgt)
msc_cs_t = np.average(msc_cs_t0[wlpkdind,:,:],axis=0,weights=wlpkdwgt)
mscb_em_t = np.average(mscb_em_t0[wlpkdind,:,:],axis=0,weights=wlpkdwgt)
mscb_cs_t = np.average(mscb_cs_t0[wlpkdind,:,:],axis=0,weights=wlpkdwgt)
mscb_em1_t = np.average(mscb_em1_t0[wlpkdind,:,:],axis=0,weights=wlpkdwgt1)
mscb_cs1_t = np.average(mscb_cs1_t0[wlpkdind,:,:],axis=0,weights=wlpkdwgt1)



mac_em  = np.zeros(volr.shape)
mac_cs  = np.zeros(volr.shape)
mscb_em  = np.zeros(volr.shape)
mscb_cs  = np.zeros(volr.shape)
mscb_em1  = np.zeros(volr.shape)
mscb_cs1  = np.zeros(volr.shape)
mec_em   = np.zeros(volr.shape)
mec_cs   = np.zeros(volr.shape)
for i,j,k in itertools.product(np.arange(volr.shape[0]),np.arange(volr.shape[1]),np.arange(volr.shape[2])):

    volpkd = volr[i,j,k] * 100
    volempkdindh = np.where( volem >= volpkd )[0][-1]
    volempkdindl = volempkdindh + 1

    mac_cs_vol2 = np.zeros([2,len(rh)])
    mscb_cs_vol2 = np.zeros([2,len(rh)])
    mscb_cs1_vol2 = np.zeros([2,len(rh)])
    msc_cs_vol2 = np.zeros([2,len(rh)])
    if ( volpkd > 98 ):
        volcspkdindl = np.where( volcs == 98 )[0][-1]
        mac_cs_vol2[0,:] = mac_cs_t[volcspkdindl,:]
        mac_cs_vol2[1,:] = mac_em_t[volempkdindh,:] # no volr=100% for
                                                             # cs, use em
        mscb_cs_vol2[0,:] = mscb_cs_t[volcspkdindl,:]
        mscb_cs_vol2[1,:] = mscb_em_t[volempkdindh,:] #
        mscb_cs1_vol2[0,:] = mscb_cs1_t[volcspkdindl,:]
        mscb_cs1_vol2[1,:] = mscb_em1_t[volempkdindh,:] 
        msc_cs_vol2[0,:] = msc_cs_t[volcspkdindl,:]
        msc_cs_vol2[1,:] = msc_em_t[volempkdindh,:] #
    else:
        volcspkdindh = np.where( volcs >= volpkd )[0][-1]
        volcspkdindl = volcspkdindh + 1

        mac_cs_vol2[0,:] = mac_cs_t[volcspkdindl,:]
        mac_cs_vol2[1,:] = mac_cs_t[volcspkdindh,:]
        mscb_cs_vol2[0,:] = mscb_cs_t[volcspkdindl,:]
        mscb_cs_vol2[1,:] = mscb_cs_t[volcspkdindh,:]
        mscb_cs1_vol2[0,:] = mscb_cs1_t[volcspkdindl,:]
        mscb_cs1_vol2[1,:] = mscb_cs1_t[volcspkdindh,:]
        msc_cs_vol2[0,:] = msc_cs_t[volcspkdindl,:]
        msc_cs_vol2[1,:] = msc_cs_t[volcspkdindh,:]


    mac_em_vol2 = mac_em_t[[volempkdindl,volempkdindh],:]
    mscb_em_vol2 = mscb_em_t[[volempkdindl,volempkdindh],:]
    mscb_em1_vol2 = mscb_em1_t[[volempkdindl,volempkdindh],:]
    msc_em_vol2 = msc_em_t[[volempkdindl,volempkdindh],:]

    volpkdl = volem[volempkdindl]
    volpkdh = volem[volempkdindh]
    volpkdr = [(volpkdh - volpkd)/(volpkdh - volpkdl), (volpkd -
        volpkdl)/(volpkdh - volpkdl)]
    #print 'volpkd = '+str(volpkd)+', volpkdl = '+str(volpkdl)+', volpkdh ='+str(volpkdh)

    mac_cs_vol = np.average(mac_cs_vol2,axis=0,weights=volpkdr)
    mac_em_vol = np.average(mac_em_vol2,axis=0,weights=volpkdr)
    mscb_cs_vol = np.average(mscb_cs_vol2,axis=0,weights=volpkdr)
    mscb_em_vol = np.average(mscb_em_vol2,axis=0,weights=volpkdr)
    mscb_cs1_vol = np.average(mscb_cs1_vol2,axis=0,weights=volpkdr)
    mscb_em1_vol = np.average(mscb_em1_vol2,axis=0,weights=volpkdr)
    msc_cs_vol = np.average(msc_cs_vol2,axis=0,weights=volpkdr)
    msc_em_vol = np.average(msc_em_vol2,axis=0,weights=volpkdr)


    rhpkd = rh_rg[i,j,k]
#    if ( rhpkd <= 30 ):
#        rhpkdindh = np.where( rh == 35 )[0][0]
#        rhpkdindl = np.where( rh == 30 )[0][0]
#        rhpkdh = 35.
#        rhpkdl = 0.
#    if ( rhpkd < 30 ):
#        rhpkdindh = np.where( rh == 30 )[0][0]
#        rhpkdindl = np.where( rh == 30 )[0][0]
#        rhpkdh = 30.
#        rhpkdl = 0.
    if ( rhpkd < 35 ):
        rhpkdindh = np.where( rh == 35 )[0][0]
        rhpkdindl = np.where( rh == 30 )[0][0]
        rhpkdh = 35.
        rhpkdl = 0.
    else:
        rhpkdindh = np.where( rh >= rhpkd )[0][0]
        rhpkdindl = rhpkdindh - 1
        rhpkdh = rh[rhpkdindh]
        rhpkdl = rh[rhpkdindl]

    mac_cs_rh = mac_cs_vol[[rhpkdindl,rhpkdindh]]
    mac_em_rh = mac_em_vol[[rhpkdindl,rhpkdindh]]
    mscb_cs_rh = mscb_cs_vol[[rhpkdindl,rhpkdindh]]
    mscb_em_rh = mscb_em_vol[[rhpkdindl,rhpkdindh]]
    mscb_cs1_rh = mscb_cs1_vol[[rhpkdindl,rhpkdindh]]
    mscb_em1_rh = mscb_em1_vol[[rhpkdindl,rhpkdindh]]
    msc_cs_rh = msc_cs_vol[[rhpkdindl,rhpkdindh]]
    msc_em_rh = msc_em_vol[[rhpkdindl,rhpkdindh]]

    rhpkdr = [(rhpkdh - rhpkd)/(rhpkdh - rhpkdl), (rhpkd -
        rhpkdl)/(rhpkdh - rhpkdl)]
    
    #print 'rhpkd = '+str(rhpkd)+', rhpkdl = '+str(rhpkdl)+', rhpkdh ='+str(rhpkdh)

    mac_cs[i,j,k] = np.average(mac_cs_rh,weights=rhpkdr)
    mac_em[i,j,k] = np.average(mac_em_rh,weights=rhpkdr)
    mscb_cs[i,j,k] = np.average(mscb_cs_rh,weights=rhpkdr)
    mscb_cs1[i,j,k] = np.average(mscb_cs1_rh,weights=rhpkdr)
    mscb_em[i,j,k] = np.average(mscb_em_rh,weights=rhpkdr)
    mscb_em1[i,j,k] = np.average(mscb_em1_rh,weights=rhpkdr)

    msc_cs = np.average(msc_cs_rh,weights=rhpkdr)
    msc_em = np.average(msc_em_rh,weights=rhpkdr)

    mec_cs[i,j,k] = mac_cs[i,j,k] + msc_cs
    mec_em[i,j,k] = mac_em[i,j,k] + msc_em

dmac = mac_cs - mac_em
dmscb = mscb_cs - mscb_em

bcso4_col = (bc_col + so4_col) * 1E3 # in units of g/m2

tau_cs = bcso4_col * mec_cs
tau_em = bcso4_col * mec_em

### One Layer Simplified Model: Calculate dRF_TOA and dRF_SFC

rs = alb_sfc/100.

a_em1  = 2. * bcso4_col * mac_em    # effective spherical absorptivity 
a_cs1  = 2. * bcso4_col * mac_cs    # considering isotropic insolation
r_em1  = 2. * bcso4_col * mscb_em1
r_cs1  = 2. * bcso4_col * mscb_cs1
t_em1  = 1. - a_em1 - r_em1
t_cs1  = 1. - a_cs1 - r_cs1
a_em  = 2. * bcso4_col * mac_em    # effective spherical absorptivity 
a_cs  = 2. * bcso4_col * mac_cs    # considering isotropic insolation
r_em  = 2. * bcso4_col * mscb_em
r_cs  = 2. * bcso4_col * mscb_cs
t_em  = 1. - a_em - r_em
t_cs  = 1. - a_cs - r_cs

# ignore MSC
#drf_sfc_smod_mon = 2.*f0 * ta * (1. - rs) * 2. * bcso4_col * dmac
#drf_toa_smod_mon = 2.*f0 * np.power(ta,2) * 2. * rs * 2. * bcso4_col * dmac

# simplified expression
#rf_sfc_em_smod_mon = -1. * f0 * ta * ((1. - rs)*a_em + \
#        r_em*np.power(1-rs,2))
#rf_sfc_cs_smod_mon = -1. * f0 * ta * ((1. - rs)*a_cs + \
#        r_cs*np.power(1-rs,2))
#rf_toa_em_smod_mon = f0 * np.power(ta,2) * (2.*rs*a_em - \
#        r_em*np.power(1-rs,2))
#rf_toa_cs_smod_mon = f0 * np.power(ta,2) * (2.*rs*a_cs - \
#        r_cs*np.power(1-rs,2))
#
#drf_sfc_smod_mon = rf_sfc_cs_smod_mon - rf_sfc_em_smod_mon
#drf_toa_smod_mon = rf_toa_cs_smod_mon - rf_toa_em_smod_mon

# full expression
rf_sfc_em_smod_mon = f0 * ta * ((1. - rs) * (t_em/(1. - rs*r_em1) - 1.))
rf_sfc_cs_smod_mon = f0 * ta * ((1. - rs) * (t_cs/(1. - rs*r_cs1) - 1.))
rf_toa_em_smod_mon = f0 * np.power(ta,2) * (rs - r_em - t_em*t_em1*rs/(1. - rs*r_em1))
rf_toa_cs_smod_mon = f0 * np.power(ta,2) * (rs - r_cs - t_cs*t_cs1*rs/(1. - rs*r_cs1))
drf_sfc_smod_mon = rf_sfc_cs_smod_mon - rf_sfc_em_smod_mon
drf_toa_smod_mon = rf_toa_cs_smod_mon - rf_toa_em_smod_mon

# Average over all months for annual mean
drf_sfc_smod = np.average(drf_sfc_smod_mon,axis=2)
drf_toa_smod = np.average(drf_toa_smod_mon,axis=2)


### Two Layer Simplified Mode: Calculate dRF_TOA and dRF_SFC
### assuming equal column density of aerosols

# Define ratio of aerosol column density split between layer 1 (bottom) and
# layer 2 (top)
acol_r_l1 = 0.5
acol_r_l2 = 0.5
# calculate t r a for layer 1
a_em_l1 = acol_r_l1 * a_em
a_cs_l1 = acol_r_l1 * a_cs
r_em_l1 = acol_r_l1 * r_em
r_cs_l1 = acol_r_l1 * r_cs
t_em_l1 = 1. - a_em_l1 - r_em_l1
t_cs_l1 = 1. - a_cs_l1 - r_cs_l1
# effective albedo of layer 1 plus surface
rs_em_l1 = r_em_l1 + np.power(t_em_l1, 2) * rs / (1. - rs * r_em_l1)
rs_cs_l1 = r_cs_l1 + np.power(t_cs_l1, 2) * rs / (1. - rs * r_cs_l1)
# calculate t r a for layer 2
a_em_l2 = acol_r_l2 * a_em
a_cs_l2 = acol_r_l2 * a_cs
r_em_l2 = acol_r_l2 * r_em
r_cs_l2 = acol_r_l2 * r_cs
t_em_l2 = 1. - a_em_l2 - r_em_l2
t_cs_l2 = 1. - a_cs_l2 - r_cs_l2

# calculate TOA forcing
rf_toa_em_smod_mon2 = f0 * np.power(ta,2) * (rs - r_em_l2 - \
        np.power(t_em_l2,2)*rs_em_l1/(1. - rs_em_l1*r_em_l2))
rf_toa_cs_smod_mon2 = f0 * np.power(ta,2) * (rs - r_cs_l2 - \
        np.power(t_cs_l2,2)*rs_cs_l1/(1. - rs_cs_l1*r_cs_l2))
drf_toa_smod_mon2 = rf_toa_cs_smod_mon2 - rf_toa_em_smod_mon2

# Average over all months for annual mean
drf_toa_smod2 = np.average(drf_toa_smod_mon2,axis=2)



### Calculate slope from historical data
for region in regions:
    
    i = regions.index(region)
    
    print region
    print 'TOA:'
    print 'EM: '
    print 'Standalone Model:  ',
    print toa_ext[i,:]
    print 'One-Layer SModel:  ',
    print np.average(rf_toa_em_smod_mon[i,:,:],axis=1)
    print 'CS: '
    print 'Standalone Model:  ',
    print toa_int[i,:]
    print 'One-Layer SModel:  ',
    print np.average(rf_toa_cs_smod_mon[i,:,:],axis=1)
    print 'dRF_TOA:  '
    print 'STandalone Model: '
    print toa_di[i,:]
    print 'One-Layer SModel:  ',
    print drf_toa_smod[i,:]
    print 'Ratio: ',
    print toa_di[i,:]/drf_toa_smod[i,:]

#    print 'Two-Layer Simplified Model:  '
#    print 'EM: ',
#    print np.average(rf_toa_em_smod_mon2,axis=2)
#    print 'CS: ',
#    print np.average(rf_toa_cs_smod_mon2,axis=2)
#    print 'dRF_TOA:  ',
#    print drf_toa_smod2[i,:]
#    print 'Ratio:             ',
#    print toa_di[i,:]/drf_toa_smod[i,:]

    print 'SFC:'
    print 'EM: '
    print 'Standalone Model:  ',
    print sfc_ext[i,:]
    print 'One-Layer SModel:  ',
    print np.average(rf_sfc_em_smod_mon[i,:,:],axis=1)
    print 'CS: '
    print 'Standalone Model:  ',
    print sfc_int[i,:]
    print 'One-Layer SModel:  ',
    print np.average(rf_sfc_cs_smod_mon[i,:,:],axis=1)
    print 'dRF_sfc:  '
    print 'STandalone Model: '
    print sfc_di[i,:]
    print 'One-Layer SModel:  ',
    print drf_sfc_smod[i,:]
    print 'Ratio:             ',
    print sfc_di[i,:]/drf_sfc_smod[i,:]
    print '\n\n'

### Write ratios and slopes to ncfile
outfile = Dataset('linreg.smodratios.nc','w',format='NETCDF4')

outfile.createDimension('region',len(regions))
outfile.createDimension('year',ratios.shape[1])
outfile.createDimension('statstype',5)

region_dim = outfile.createVariable('region',np.int32,('region',))
year_dim   = outfile.createVariable('year',np.int32,('year',))
stats_dim  = outfile.createVariable('statstype',np.int32,('statstype'))

description = ''
for region, rgid in zip(regions, rgids):
    description += region + ': ' + str(rgid) + ',  '
region_dim.description =description

description = '0-slope; 1-intercept; 2-r_value; 3-p_value; 4-std_err'
stats_dim.description =description

stats_out = outfile.createVariable('regstats',np.float32,('region','statstype'))

ratios_out = outfile.createVariable('smodratios',np.float32,('region','year'))

region_dim[:] = rgids
year_dim[:]   = years[:-1]
stats_dim[:]  = [0,1,2,3,4]

stats_out[:]  = regstats
ratios_out[:] = ratios

outfile.close()

### Calculate correlation between slopes and ratios
#slope, intercept, r_value, p_value, std_err = stats.linregress(np.append(ratios,[0]),np.append(slopes,[0]))
#slope, intercept, r_value, p_value, std_err = stats.linregress(ratios,slopes)
#
#print 'correlation between slopes and ratios:'
#print 'slope = ' + str(slope),
#print 'intercept = ' + str(intercept),
#print 'r_square = ' + str(np.power(r_value,2))

#plt.scatter(slopes, ratios, c='red', cmap='RdYlBu', s=150, alpha=0.8,linewidths=0.5, marker='o')
#plt.grid()
#plt.show()
