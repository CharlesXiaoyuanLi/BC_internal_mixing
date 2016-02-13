### Calculate the linear regression of TOA vs. SFC

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

toa_di = fh_toa.variables['rf_di'][:][:][:]
sfc_di = fh_sfc.variables['rf_di'][:][:][:]
f0      = fh_par.variables['f0'][:][:][:]
ta      = fh_par.variables['ta'][:][:][:]
alb_sfc = fh_par.variables['alb_sfc'][:][:][:]

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

bc_col = np.average(bc_col_t[:,:-1,:],axis=2)
so4_col = np.average(so4_col_t[:,:-1,:],axis=2)

volr = so4_col/1.74/(so4_col/1.74 + bc_col/1.8)

print volr.shape

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

mac_em_t = fh_em.variables['beta_e'][:][:][:] * (1. - \
        np.array(fh_em.variables['ssa'][:][:][:]))
mac_cs_t = fh_cs.variables['beta_e'][:][:][:] * (1. - \
        np.array(fh_cs.variables['ssa'][:][:][:]))

msc_em_t = fh_em.variables['beta_e'][:][:][:] * fh_em.variables['ssa'][:][:][:]
g_em   = fh_em.variables['g'][:][:][:]
msc_cs_t = fh_cs.variables['beta_e'][:][:][:] * fh_cs.variables['ssa'][:][:][:]
g_cs   = fh_cs.variables['g'][:][:][:]

mscb_em_t = msc_em_t * 0.5 * (1. - g_em)  # currently not appplied: effective backscattering ratio
mscb_cs_t = msc_cs_t * 0.5 * (1. - g_cs)  # by considering isotropic insolation

# Pick the MAC data that match
wlpkd = 0.55
wlpkdind = np.where((wl >= wlpkd-0.02) & (wl <= wlpkd+0.02))[0].item()

rhpkd = 70.
rhpkdind = np.where( rh == rhpkd )[0].item()

mac_em  = np.zeros(volr.shape)
mac_cs  = np.zeros(volr.shape)
mscb_em  = np.zeros(volr.shape)
mscb_cs  = np.zeros(volr.shape)
for i,j in itertools.product(np.arange(volr.shape[0]),np.arange(volr.shape[1])):
    volpkd = volr[i,j] * 100
    if ( volpkd >= 80 and volpkd <=99):
        volempkdind = np.where( (volem < volpkd + 1) & (volem >= volpkd - \
            1))[0].item()
        volcspkdind = np.where( (volcs < volpkd + 1) & (volcs >= volpkd - \
                1))[0].item()
    elif ( volpkd > 99 ):
        volcspkdind = np.where( volcs == 98 )
        volempkdind = np.where( volem == 98 )
    elif ( volpkd > 79.5 and volpkd < 80 ):
        volempkdind = np.where( volem == 80 )
        volcspkdind = np.where( volcs == 80 )
    else:
        volempkdind = np.where( (volem < volpkd + 2.5) & (volem >= volpkd - \
            2.5))[0].item()
        volcspkdind = np.where( (volcs < volpkd + 2.5) & (volcs >= volpkd - \
            2.5))[0].item()
    mac_cs[i,j] = mac_cs_t[wlpkdind,volcspkdind,rhpkdind]
    mac_em[i,j] = mac_em_t[wlpkdind,volempkdind,rhpkdind]
    mscb_cs[i,j] = mscb_cs_t[wlpkdind,volcspkdind,rhpkdind]
    mscb_em[i,j] = mscb_em_t[wlpkdind,volempkdind,rhpkdind]

dmac = mac_cs - mac_em
dmscb = mscb_cs - mscb_em

bcso4_col = (bc_col + so4_col) * 1E3 # in units of g/m2


### One Layer Simplified Model: Calculate dRF_TOA and dRF_SFC

rs = np.average(alb_sfc,axis=2)/100.

a_em  = 2. * bcso4_col * mac_em    # effective spherical absorptivity 
a_cs  = 2. * bcso4_col * mac_cs    # considering isotropic insolation
r_em  = 2. * bcso4_col * mscb_em
r_cs  = 2. * bcso4_col * mscb_cs
t_em  = 1. - a_em - r_em
t_cs  = 1. - a_cs - r_cs

# ignore MSC
#drf_sfc_smod_mon = f0 * ta * (1. - alb_sfc/100.) * bcso4_col * dmac
#drf_toa_smod_mon = f0 * np.power(ta,2) * 2. * alb_sfc/100. * bcso4_col * dmac

# simplified expression
drf_sfc_smod = -1. * np.average(f0,axis=2) * np.average(ta,axis=2) * ((1. - rs)*2*dmac*bcso4_col + \
        2*dmscb*bcso4_col*np.power(1-rs,2)) 
drf_toa_smod = np.average(f0,axis=2) * np.power(np.average(ta,axis=2),2) * (4.*rs*dmac*bcso4_col - \
        2*dmscb*bcso4_col*np.power(1-rs,2))


# full expression
#rf_sfc_em_smod_mon = f0 * ta * ((1. - rs) * (t_em/(1. - rs*r_em) - 1.))
#rf_sfc_cs_smod_mon = f0 * ta * ((1. - rs) * (t_cs/(1. - rs*r_cs) - 1.))
#rf_toa_em_smod_mon = f0 * np.power(ta,2) * (rs - r_em - np.power(t_em,2)*rs/(1. - rs*r_em))
#rf_toa_cs_smod_mon = f0 * np.power(ta,2) * (rs - r_cs - np.power(t_cs,2)*rs/(1. - rs*r_cs))
#drf_sfc_smod_mon = rf_sfc_cs_smod_mon - rf_sfc_em_smod_mon
#drf_toa_smod_mon = rf_toa_cs_smod_mon - rf_toa_em_smod_mon


### Two Layer Simplified Mode: Calculate dRF_TOA and dRF_SFC
### assuming equal column density of aerosols

# Define ratio of aerosol column density split between layer 1 (bottom) and
# layer 2 (top)
acol_r_l1 = 0.8
acol_r_l2 = 0.2
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
rf_toa_em_smod_mon2 = np.average(f0,axis=2) * np.power(np.average(ta,axis=2),2) * (rs - r_em_l2 - \
        np.power(t_em_l2,2)*rs_em_l1/(1. - rs_em_l1*r_em_l2))
rf_toa_cs_smod_mon2 = np.average(f0,axis=2) * np.power(np.average(ta,axis=2),2) * (rs - r_cs_l2 - \
        np.power(t_cs_l2,2)*rs_cs_l1/(1. - rs_cs_l1*r_cs_l2))
drf_toa_smod2 = rf_toa_cs_smod_mon2 - rf_toa_em_smod_mon2




### Calculate slope from historical data
for region in regions:
    
    i = regions.index(region)
    
    print region
    print 'TOA:'
    print 'Standalone Model:  ',
    print toa_di[i,:]
    print 'One-Layer Simplified Model:  '
#    print 'EM: ',
#    print np.average(rf_toa_em_smod_mon,axis=2)
#    print 'CS: ',
#    print np.average(rf_toa_cs_smod_mon,axis=2)
#    print 'dRF_TOA:  ',
    print drf_toa_smod[i,:]
    print 'Two-Layer Simplified Model:  '
#    print 'EM: ',
#    print np.average(rf_toa_em_smod_mon2,axis=2)
#    print 'CS: ',
#    print np.average(rf_toa_cs_smod_mon2,axis=2)
#    print 'dRF_TOA:  ',
    print drf_toa_smod2[i,:]
    print 'Ratio:             ',
    print toa_di[i,:]/drf_toa_smod[i,:]
    print 'SFC:'
    print 'Standalone Model:  ',
    print sfc_di[i,:]
    print 'Simplified Model:  ',
    print drf_sfc_smod[i,:]
    print 'Ratio:             ',
    print sfc_di[i,:]/drf_sfc_smod[i,:]
    print '\n\n'

sys.exit()

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
