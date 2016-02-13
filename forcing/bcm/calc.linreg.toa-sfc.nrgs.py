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


gr = 1.
if (imtype == 'imbcm_ndist'):
    drtda = -0.2
elif (imtype == 'csbcm_ndist'):
    drtda = 0


### Read data

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

bc_col = bc_col_t[:,:-1,:]
so4_col = so4_col_t[:,:-1,:]

volr = so4_col/1.74/(so4_col/1.74 + bc_col/1.8)

#### Calculate average weight of bc and so4 col den for each month as an example
#bcso4_col_t = bc_col + so4_col
#
## convert to weight for each year
#bcso4_col_wgt_t = np.array([ bcso4_col_t[:,:,i] / np.sum(bcso4_col_t,axis=2) for i in
#        np.arange(12) ])  # shape = (12,15,8)
## average the weight over all years for each month
#bcso4_col_avg_wgt_t = np.average( bcso4_col_wgt_t, axis=2 ) # shape = (12,15)
#bcso4_col_avg_wgt = np.transpose( bcso4_col_avg_wgt_t )
#
#print bcso4_col_avg_wgt.shape

### Read in MAC & MSC data
ncfile_em = '../../miedata/sul-bcm-em-0.0118.nc'
if (imtype == 'imbcm_ndist'):
    ncfile_cs = '../../miedata/sul-bcm-im-ndist.nc'
elif (imtype == 'csbcm_ndist'):
    ncfile_cs = '../../miedata/sul-bcm-coat-ndist.new.nc'

fh_em = Dataset(ncfile_em, 'r')
fh_cs = Dataset(ncfile_cs, 'r')

wl      = fh_em.variables['wl'][:]
volem   = fh_em.variables['vol'][:]
rh      = fh_em.variables['RH'][:]
volcs   = fh_cs.variables['vol'][:]

mac_em = fh_em.variables['beta_e'][:][:][:] * (1. - \
        np.array(fh_em.variables['ssa'][:][:][:]))
mac_cs = fh_cs.variables['beta_e'][:][:][:] * (1. - \
        np.array(fh_cs.variables['ssa'][:][:][:]))

msc_em = fh_em.variables['beta_e'][:][:][:] * fh_em.variables['ssa'][:][:][:]
g_em   = fh_em.variables['g'][:][:][:]
msc_cs = fh_cs.variables['beta_e'][:][:][:] * fh_cs.variables['ssa'][:][:][:]
g_cs   = fh_cs.variables['g'][:][:][:]

mscb_em = msc_em * 0.5 * (1. - gr*g_em)
mscb_cs = msc_cs * 0.5 * (1. - gr*g_cs)

# Pick the MAC data that match
wlpkd = 0.55
wlpkdind = np.where((wl >= wlpkd-0.02) & (wl <= wlpkd+0.02))[0].item()

rhpkd = 60.
rhpkdind = np.where( rh == rhpkd )[0].item()

dmac  = np.zeros(volr.shape)
dmscb = np.zeros(volr.shape)
for i,j,k in itertools.product(np.arange(volr.shape[0]),np.arange(volr.shape[1]),np.arange(volr.shape[2])):
    volpkd = volr[i,j,k] * 100
    if ( volpkd >= 80 and volpkd <=99):
        volempkdind = np.where( (volem < volpkd + 1) & (volem >= volpkd - \
            1))[0].item()
        volcspkdind = np.where( (volcs < volpkd + 1) & (volcs >= volpkd - \
                1))[0].item()
    elif ( volpkd > 99 ):
        volcspkdind = np.where( volcs == 98 )
        volempkdind = np.where( volem == 98 )
    elif ( volpkd >79.5 and volpkd < 80):
        volempkdind = np.where( volem == 80 )
        volcspkdind = np.where( volcs == 80 )
    else:
        volempkdind = np.where( (volem < volpkd + 2.5) & (volem >= volpkd - \
            2.5))[0].item()
        volcspkdind = np.where( (volcs < volpkd + 2.5) & (volcs >= volpkd - \
            2.5))[0].item()
    dmac[i,j,k] = mac_cs[wlpkdind,volcspkdind,rhpkdind] - \
        mac_em[wlpkdind,volempkdind,rhpkdind]
    dmscb[i,j,k] = mscb_cs[wlpkdind,volcspkdind,rhpkdind] - \
        mscb_em[wlpkdind,volempkdind,rhpkdind]


bcso4_col = (bc_col + so4_col) * 1E3 # in units of g/m2

alpha = np.zeros(bcso4_col.shape)
#avgang = np.array([45,40,40,35,45,45,40,45,45,40,40,50,40,40,65])
avgang = np.array([45,45,30,35,45,40,30,55,45,40,35,45,30,35,70])
for i in np.arange(alpha.shape[0]):
    alpha[i,:,:] = 1./np.cos(avgang[i]*1./180*np.pi)


###  Calculate ratio from atmospheric parameters
###     according to the simplified model

#fact = 1.
# Ignore MAC and MSC since their seasonal variability is negligible compared
# with that of column density
#fact = bcso4_col 
#nume = f0 * ta * (1. - alb_sfc/100.) * fact
#deno = f0 * np.power(ta,2) * 2. * alb_sfc/100. * fact

#nume = f0 * ta * ((1. - alb_sfc/100.) - drtda*np.power(1-alb_sfc/100.,2)) * fact
#deno = f0 * np.power(ta,2) * (2. * alb_sfc/100. + drtda*np.power(1-alb_sfc/100.,2)) * fact

# ignore only MSC
#fact = bcso4_col * dmac
#nume = f0 * ta * (1. - alb_sfc/100.) * fact
#deno = f0 * np.power(ta,2) * 2. * alb_sfc/100. * fact

# full expression
#nume = f0 * ta * ((1. - alb_sfc/100.)*dmac*bcso4_col + \
#                dmscb*bcso4_col*np.power(1-alb_sfc/100.,2)) 
#deno = f0 * np.power(ta,2) * (2.*alb_sfc/100.*dmac*bcso4_col - \
#                dmscb*bcso4_col*np.power(1-alb_sfc/100.,2))

# improved simplified expression
#fact = bcso4_col
#rs = alb_sfc/100.
#nume = f0 * ta * ((1.-rs)*dmac*alpha + drtda * dmscb*(1.-rs)*(alpha-rs)) * fact
#deno = f0 * ta*ta * (rs*dmac*(1.+alpha) - drtda * dmscb*(rs*rs + alpha - rs*(1.+alpha))) * fact

# improved simplified expression without MSC
fact = bcso4_col
rs = alb_sfc/100.
nume = f0 * ta * ((1.-rs)*alpha + drtda * (1.-rs)*(alpha-rs)) * fact
deno = f0 * ta*ta * (rs*(1.+alpha) - drtda * (rs*rs + alpha - rs*(1.+alpha))) * fact

ratios = np.average(nume, axis=2)/np.average(deno, axis=2)

### Calculate slope from historical data
regstats = np.zeros((ratios.shape[0],5)) # 5 variables for stats; shape of ratios = (15, nyr)
for region in regions:
    
    i = regions.index(region)
    slope, intercept, r_value, p_value, std_err = stats.linregress(toa_di[i,:],
            sfc_di[i,:]*-1.)
    regstats[i,0] = slope
    regstats[i,1] = intercept
    regstats[i,2] = r_value
    regstats[i,3] = p_value
    regstats[i,4] = std_err
    
    print region + ': '
    print 'slope = ' + str(slope),
    print 'intercept = ' + str(intercept),
    print 'r_square = ' + str(np.power(r_value,2))

    print 'Derived Ratio = '
    for ra in ratios[i,:]:
        print str(ra)+'\t',

    print '\n\n'

### Write ratios and slopes to ncfile
outfile = Dataset('linreg.smodratios.'+imtype+'.nc','w',format='NETCDF4')

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
