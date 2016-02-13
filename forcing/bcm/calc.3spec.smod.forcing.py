### Calculate the forcing according to Simplified RTM for 3-species mixing
#(BC, sulfate and OC)

import numpy as np
from netCDF4 import Dataset
from scipy import stats
import itertools
import matplotlib.pyplot as plt
import sys

### OM/OC ratio

omr = 1.7

### Read bcsulocem from Standalone Model

fn_toa = "toa_embcm_imbcm_ndist_15rgs"+".nc"
fn_sfc = "sfc_embcm_imbcm_ndist_15rgs"+".nc"
fn_par = "atmos_param_15rgs"+".nc"

fh_toa = Dataset(fn_toa, 'r')
fh_sfc = Dataset(fn_sfc, 'r')
fh_par = Dataset(fn_par, 'r')

toa_bcsuloc_em = fh_toa.variables['rf_bcsuloc_em'][:][:]
sfc_bcsuloc_em = fh_sfc.variables['rf_bcsuloc_em'][:][:]
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
fn_col = "./aerocol_15rgs.nc"
fh_col = Dataset(fn_col, 'r')

bc_col_t  = fh_col.variables['bc_col'][:][:][:] #no 2002 
so4_col_t = fh_col.variables['so4_col'][:][:][:] 
oc_col_t  = fh_col.variables['oc_col'][:][:][:]

years = fh_col.variables['year'][:]
rgids = fh_col.variables['region'][:]

bc_col = bc_col_t[:,:-1,:]
so4_col = so4_col_t[:,:-1,:]
oc_col = omr*oc_col_t[:,:-1,:]

vols = so4_col/1.74/(so4_col/1.74 + bc_col/1.8)
volo = oc_col/1.5/(oc_col/1.5 + bc_col/1.8)

### Read in MAC & MSC data
ncfile_3em = '../../miedata/bcm-sul-oc-em.nc'
ncfile_2im = '../../miedata/bcm-sul-im-oc-em.nc'
ncfile_3im = '../../miedata/bcm-sul-oc-im.nc'
ncfile_2cs = '../../miedata/bcm-sul-coat-oc-em.nc'
ncfile_3cs = '../../miedata/bcm-sul-oc-coat.nc'

fh_3em = Dataset(ncfile_3em, 'r')
fh_2im = Dataset(ncfile_2im, 'r')
fh_3im = Dataset(ncfile_3im, 'r')
fh_2cs = Dataset(ncfile_2cs, 'r')
fh_3cs = Dataset(ncfile_3cs, 'r')

wl      = fh_3em.variables['wl'][:]
vol     = fh_3em.variables['vols'][:]
rh      = fh_3em.variables['RH'][:]

mac_3em_t0 = fh_3em.variables['beta_e'][:][:][:][:] * (1. - \
        np.array(fh_3em.variables['ssa'][:][:][:][:]))
mac_2im_t0 = fh_2im.variables['beta_e'][:][:][:][:] * (1. - \
        np.array(fh_2im.variables['ssa'][:][:][:][:]))
mac_3im_t0 = fh_3im.variables['beta_e'][:][:][:][:] * (1. - \
        np.array(fh_3im.variables['ssa'][:][:][:][:]))
mac_2cs_t0 = fh_2cs.variables['beta_e'][:][:][:][:] * (1. - \
        np.array(fh_2cs.variables['ssa'][:][:][:][:]))
mac_3cs_t0 = fh_3cs.variables['beta_e'][:][:][:][:] * (1. - \
        np.array(fh_3cs.variables['ssa'][:][:][:][:]))

msc_3em_t0 = fh_3em.variables['beta_e'][:][:][:][:] * fh_3em.variables['ssa'][:][:][:][:]
g_3em   = fh_3em.variables['g'][:][:][:][:]
msc_2im_t0 = fh_2im.variables['beta_e'][:][:][:][:] * fh_2im.variables['ssa'][:][:][:][:]
g_2im   = fh_2im.variables['g'][:][:][:][:]
msc_3im_t0 = fh_3im.variables['beta_e'][:][:][:][:] * fh_3im.variables['ssa'][:][:][:][:]
g_3im   = fh_3im.variables['g'][:][:][:][:]
msc_2cs_t0 = fh_2cs.variables['beta_e'][:][:][:][:] * fh_2cs.variables['ssa'][:][:][:][:]
g_2cs   = fh_2cs.variables['g'][:][:][:][:]
msc_3cs_t0 = fh_3cs.variables['beta_e'][:][:][:][:] * fh_3cs.variables['ssa'][:][:][:][:]
g_3cs   = fh_3cs.variables['g'][:][:][:][:]

mscb_3em_t0 = msc_3em_t0 * 0.5 * (1. - g_3em)
mscb_2im_t0 = msc_2im_t0 * 0.5 * (1. - g_2im)
mscb_3im_t0 = msc_3im_t0 * 0.5 * (1. - g_3im)
mscb_2cs_t0 = msc_2cs_t0 * 0.5 * (1. - g_2cs)
mscb_3cs_t0 = msc_3cs_t0 * 0.5 * (1. - g_3cs)

fh_3em.close() 
fh_2im.close() 
fh_3im.close() 
fh_2cs.close() 
fh_3cs.close() 


# Pick the MAC data that match
wlpkdind = np.where( wl < 2.5 )[0]
wlpkdwgt = np.array([0.1*250,0.2*375,0.3*500,0.8*430,1.6*255,1.9*180,1.4*180,0.7*80,0])
#wlpkdwgt = np.array([0,0,0,0,0,1,0,0,0])
#wlpkdwgt = np.array([0*250,0*375,0.25*500,0.45*430,0.8*255,1.25*180,0.3*180,0*80,0])

print wlpkdind
print wl[wlpkdind]

mac_3em_t = np.average(mac_3em_t0[wlpkdind,:,:,:],axis=0,weights=wlpkdwgt)
msc_3em_t = np.average(msc_3em_t0[wlpkdind,:,:,:],axis=0,weights=wlpkdwgt)
mscb_3em_t = np.average(mscb_3em_t0[wlpkdind,:,:,:],axis=0,weights=wlpkdwgt)
mac_2im_t = np.average(mac_2im_t0[wlpkdind,:,:,:],axis=0,weights=wlpkdwgt)
msc_2im_t = np.average(msc_2im_t0[wlpkdind,:,:,:],axis=0,weights=wlpkdwgt)
mscb_2im_t = np.average(mscb_2im_t0[wlpkdind,:,:,:],axis=0,weights=wlpkdwgt)
mac_3im_t = np.average(mac_3im_t0[wlpkdind,:,:,:],axis=0,weights=wlpkdwgt)
msc_3im_t = np.average(msc_3im_t0[wlpkdind,:,:,:],axis=0,weights=wlpkdwgt)
mscb_3im_t = np.average(mscb_3im_t0[wlpkdind,:,:,:],axis=0,weights=wlpkdwgt)
mac_2cs_t = np.average(mac_2cs_t0[wlpkdind,:,:,:],axis=0,weights=wlpkdwgt)
msc_2cs_t = np.average(msc_2cs_t0[wlpkdind,:,:,:],axis=0,weights=wlpkdwgt)
mscb_2cs_t = np.average(mscb_2cs_t0[wlpkdind,:,:,:],axis=0,weights=wlpkdwgt)
mac_3cs_t = np.average(mac_3cs_t0[wlpkdind,:,:,:],axis=0,weights=wlpkdwgt)
msc_3cs_t = np.average(msc_3cs_t0[wlpkdind,:,:,:],axis=0,weights=wlpkdwgt)
mscb_3cs_t = np.average(mscb_3cs_t0[wlpkdind,:,:,:],axis=0,weights=wlpkdwgt)


mac_3em  = np.zeros(vols.shape)
mac_2im  = np.zeros(vols.shape)
mac_3im  = np.zeros(vols.shape)
mac_2cs  = np.zeros(vols.shape)
mac_3cs  = np.zeros(vols.shape)

mscb_3em  = np.zeros(vols.shape)
mscb_2im  = np.zeros(vols.shape)
mscb_3im  = np.zeros(vols.shape)
mscb_2cs  = np.zeros(vols.shape)
mscb_3cs  = np.zeros(vols.shape)

mec_3em   = np.zeros(vols.shape)
mec_2im   = np.zeros(vols.shape)
mec_3im   = np.zeros(vols.shape)
mec_2cs   = np.zeros(vols.shape)
mec_3cs   = np.zeros(vols.shape)

for i,j,k in itertools.product(np.arange(vols.shape[0]),np.arange(vols.shape[1]),np.arange(vols.shape[2])):

    volspkd = vols[i,j,k] * 100

    mac_3em_vol2 = np.zeros([2,len(rh)])
    mscb_3em_vol2 = np.zeros([2,len(rh)])
    msc_3em_vol2 = np.zeros([2,len(rh)])

    mac_2im_vol2 = np.zeros([2,len(rh)])
    mscb_2im_vol2 = np.zeros([2,len(rh)])
    msc_2im_vol2 = np.zeros([2,len(rh)])

    mac_3im_vol2 = np.zeros([2,len(rh)])
    mscb_3im_vol2 = np.zeros([2,len(rh)])
    msc_3im_vol2 = np.zeros([2,len(rh)])

    mac_2cs_vol2 = np.zeros([2,len(rh)])
    mscb_2cs_vol2 = np.zeros([2,len(rh)])
    msc_2cs_vol2 = np.zeros([2,len(rh)])

    mac_3cs_vol2 = np.zeros([2,len(rh)])
    mscb_3cs_vol2 = np.zeros([2,len(rh)])
    msc_3cs_vol2 = np.zeros([2,len(rh)])

    if ( volspkd >= 98 ):
        volspkdindh = np.where( vol == 98 )[0][0]
        volspkdindl = volspkdindh
    else:
        volspkdindh = np.where( vol >= volspkd )[0][-1]
        volspkdindl = volspkdindh + 1

    mac_3em_vol2 = mac_3em_t[[volspkdindl,volspkdindh],:,:]
    mscb_3em_vol2 = mscb_3em_t[[volspkdindl,volspkdindh],:,:]
    msc_3em_vol2 = msc_3em_t[[volspkdindl,volspkdindh],:,:]

    mac_2im_vol2 = mac_2im_t[[volspkdindl,volspkdindh],:,:]
    mscb_2im_vol2 = mscb_2im_t[[volspkdindl,volspkdindh],:,:]
    msc_2im_vol2 = msc_2im_t[[volspkdindl,volspkdindh],:,:]

    mac_3im_vol2 = mac_3im_t[[volspkdindl,volspkdindh],:,:]
    mscb_3im_vol2 = mscb_3im_t[[volspkdindl,volspkdindh],:,:]
    msc_3im_vol2 = msc_3im_t[[volspkdindl,volspkdindh],:,:]

    mac_2cs_vol2 = mac_2cs_t[[volspkdindl,volspkdindh],:,:]
    mscb_2cs_vol2 = mscb_2cs_t[[volspkdindl,volspkdindh],:,:]
    msc_2cs_vol2 = msc_2cs_t[[volspkdindl,volspkdindh],:,:]

    mac_3cs_vol2 = mac_3cs_t[[volspkdindl,volspkdindh],:,:]
    mscb_3cs_vol2 = mscb_3cs_t[[volspkdindl,volspkdindh],:,:]
    msc_3cs_vol2 = msc_3cs_t[[volspkdindl,volspkdindh],:,:]

    if ( vol[volspkdindl] == 98 ):
        volspkdr = [1,1]
    else:
        volspkdl = vol[volspkdindl]
        volspkdh = vol[volspkdindh]
        volspkdr = [(volspkdh - volspkd)/(volspkdh - volspkdl), (volspkd -
            volspkdl)/(volspkdh - volspkdl)]

    #print 'volspkd = '+str(volspkd)+', volspkdl = '+str(volspkdl)+', volspkdh ='+str(volspkdh)

    mac_3em_vols = np.average(mac_3em_vol2,axis=0,weights=volspkdr)
    mscb_3em_vols = np.average(mscb_3em_vol2,axis=0,weights=volspkdr)
    msc_3em_vols = np.average(msc_3em_vol2,axis=0,weights=volspkdr)

    mac_2im_vols = np.average(mac_2im_vol2,axis=0,weights=volspkdr)
    mscb_2im_vols = np.average(mscb_2im_vol2,axis=0,weights=volspkdr)
    msc_2im_vols = np.average(msc_2im_vol2,axis=0,weights=volspkdr)

    mac_3im_vols = np.average(mac_3im_vol2,axis=0,weights=volspkdr)
    mscb_3im_vols = np.average(mscb_3im_vol2,axis=0,weights=volspkdr)
    msc_3im_vols = np.average(msc_3im_vol2,axis=0,weights=volspkdr)

    mac_2cs_vols = np.average(mac_2cs_vol2,axis=0,weights=volspkdr)
    mscb_2cs_vols = np.average(mscb_2cs_vol2,axis=0,weights=volspkdr)
    msc_2cs_vols = np.average(msc_2cs_vol2,axis=0,weights=volspkdr)

    mac_3cs_vols = np.average(mac_3cs_vol2,axis=0,weights=volspkdr)
    mscb_3cs_vols = np.average(mscb_3cs_vol2,axis=0,weights=volspkdr)
    msc_3cs_vols = np.average(msc_3cs_vol2,axis=0,weights=volspkdr)



    volopkd = volo[i,j,k] * 100

    if ( volopkd >= 98 ):
        volopkdindh = np.where( vol == 98 )[0][0]
        volopkdindl = volopkdindh
    else:
        volopkdindh = np.where( vol >= volopkd )[0][-1]
        volopkdindl = volopkdindh + 1

    mac_3em_vol2 = mac_3em_vols[[volopkdindl,volopkdindh],:]
    mscb_3em_vol2 = mscb_3em_vols[[volopkdindl,volopkdindh],:]
    msc_3em_vol2 = msc_3em_vols[[volopkdindl,volopkdindh],:]

    mac_2im_vol2 = mac_2im_vols[[volopkdindl,volopkdindh],:]
    mscb_2im_vol2 = mscb_2im_vols[[volopkdindl,volopkdindh],:]
    msc_2im_vol2 = msc_2im_vols[[volopkdindl,volopkdindh],:]

    mac_3im_vol2 = mac_3im_vols[[volopkdindl,volopkdindh],:]
    mscb_3im_vol2 = mscb_3im_vols[[volopkdindl,volopkdindh],:]
    msc_3im_vol2 = msc_3im_vols[[volopkdindl,volopkdindh],:]

    mac_2cs_vol2 = mac_2cs_vols[[volopkdindl,volopkdindh],:]
    mscb_2cs_vol2 = mscb_2cs_vols[[volopkdindl,volopkdindh],:]
    msc_2cs_vol2 = msc_2cs_vols[[volopkdindl,volopkdindh],:]

    mac_3cs_vol2 = mac_3cs_vols[[volopkdindl,volopkdindh],:]
    mscb_3cs_vol2 = mscb_3cs_vols[[volopkdindl,volopkdindh],:]
    msc_3cs_vol2 = msc_3cs_vols[[volopkdindl,volopkdindh],:]

    if ( vol[volopkdindl] == 98 ):
        volopkdr = [1,1]
    else:
        volopkdl = vol[volopkdindl]
        volopkdh = vol[volopkdindh]
        volopkdr = [(volopkdh - volopkd)/(volopkdh - volopkdl), (volopkd -
            volopkdl)/(volopkdh - volopkdl)]

    mac_3em_vol = np.average(mac_3em_vol2,axis=0,weights=volopkdr)
    mscb_3em_vol = np.average(mscb_3em_vol2,axis=0,weights=volopkdr)
    msc_3em_vol = np.average(msc_3em_vol2,axis=0,weights=volopkdr)

    mac_2im_vol = np.average(mac_2im_vol2,axis=0,weights=volopkdr)
    mscb_2im_vol = np.average(mscb_2im_vol2,axis=0,weights=volopkdr)
    msc_2im_vol = np.average(msc_2im_vol2,axis=0,weights=volopkdr)

    mac_3im_vol = np.average(mac_3im_vol2,axis=0,weights=volopkdr)
    mscb_3im_vol = np.average(mscb_3im_vol2,axis=0,weights=volopkdr)
    msc_3im_vol = np.average(msc_3im_vol2,axis=0,weights=volopkdr)

    mac_2cs_vol = np.average(mac_2cs_vol2,axis=0,weights=volopkdr)
    mscb_2cs_vol = np.average(mscb_2cs_vol2,axis=0,weights=volopkdr)
    msc_2cs_vol = np.average(msc_2cs_vol2,axis=0,weights=volopkdr)

    mac_3cs_vol = np.average(mac_3cs_vol2,axis=0,weights=volopkdr)
    mscb_3cs_vol = np.average(mscb_3cs_vol2,axis=0,weights=volopkdr)
    msc_3cs_vol = np.average(msc_3cs_vol2,axis=0,weights=volopkdr)


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

    mac_3em_rh = mac_3em_vol[[rhpkdindl,rhpkdindh]]
    mscb_3em_rh = mscb_3em_vol[[rhpkdindl,rhpkdindh]]
    msc_3em_rh = msc_3em_vol[[rhpkdindl,rhpkdindh]]

    mac_2im_rh = mac_2im_vol[[rhpkdindl,rhpkdindh]]
    mscb_2im_rh = mscb_2im_vol[[rhpkdindl,rhpkdindh]]
    msc_2im_rh = msc_2im_vol[[rhpkdindl,rhpkdindh]]

    mac_3im_rh = mac_3im_vol[[rhpkdindl,rhpkdindh]]
    mscb_3im_rh = mscb_3im_vol[[rhpkdindl,rhpkdindh]]
    msc_3im_rh = msc_3im_vol[[rhpkdindl,rhpkdindh]]

    mac_2cs_rh = mac_2cs_vol[[rhpkdindl,rhpkdindh]]
    mscb_2cs_rh = mscb_2cs_vol[[rhpkdindl,rhpkdindh]]
    msc_2cs_rh = msc_2cs_vol[[rhpkdindl,rhpkdindh]]

    mac_3cs_rh = mac_3cs_vol[[rhpkdindl,rhpkdindh]]
    mscb_3cs_rh = mscb_3cs_vol[[rhpkdindl,rhpkdindh]]
    msc_3cs_rh = msc_3cs_vol[[rhpkdindl,rhpkdindh]]

    rhpkdr = [(rhpkdh - rhpkd)/(rhpkdh - rhpkdl), (rhpkd -
        rhpkdl)/(rhpkdh - rhpkdl)]
    
    mac_3em[i,j,k] = np.average(mac_3em_rh,weights=rhpkdr)
    mscb_3em[i,j,k] = np.average(mscb_3em_rh,weights=rhpkdr)
    msc_3em = np.average(msc_3em_rh,weights=rhpkdr)
    mec_3em[i,j,k] = mac_3em[i,j,k] + msc_3em

    mac_2im[i,j,k] = np.average(mac_2im_rh,weights=rhpkdr)
    mscb_2im[i,j,k] = np.average(mscb_2im_rh,weights=rhpkdr)
    msc_2im = np.average(msc_2im_rh,weights=rhpkdr)
    mec_2im[i,j,k] = mac_2im[i,j,k] + msc_2im

    mac_3im[i,j,k] = np.average(mac_3im_rh,weights=rhpkdr)
    mscb_3im[i,j,k] = np.average(mscb_3im_rh,weights=rhpkdr)
    msc_3im = np.average(msc_3im_rh,weights=rhpkdr)
    mec_3im[i,j,k] = mac_3im[i,j,k] + msc_3im

    mac_2cs[i,j,k] = np.average(mac_2cs_rh,weights=rhpkdr)
    mscb_2cs[i,j,k] = np.average(mscb_2cs_rh,weights=rhpkdr)
    msc_2cs = np.average(msc_2cs_rh,weights=rhpkdr)
    mec_2cs[i,j,k] = mac_2cs[i,j,k] + msc_2cs

    mac_3cs[i,j,k] = np.average(mac_3cs_rh,weights=rhpkdr)
    mscb_3cs[i,j,k] = np.average(mscb_3cs_rh,weights=rhpkdr)
    msc_3cs = np.average(msc_3cs_rh,weights=rhpkdr)
    mec_3cs[i,j,k] = mac_3cs[i,j,k] + msc_3cs

bcso4oc_col = (bc_col + so4_col + oc_col) * 1E3 # in units of g/m2

tau_3em = bcso4oc_col * mec_3em
tau_2im = bcso4oc_col * mec_2im
tau_3im = bcso4oc_col * mec_3im
tau_2cs = bcso4oc_col * mec_2cs
tau_3cs = bcso4oc_col * mec_3cs


### One Layer Simplified Model: Calculate dRF_TOA and dRF_SFC

rs = alb_sfc/100.

a_3em  = 2. * bcso4oc_col * mac_3em    # effective spherical absorptivity 
r_3em  = 2. * bcso4oc_col * mscb_3em
t_3em  = 1. - a_3em - r_3em

a_2im  = 2. * bcso4oc_col * mac_2im    # effective spherical absorptivity 
r_2im  = 2. * bcso4oc_col * mscb_2im
t_2im  = 1. - a_2im - r_2im

a_3im  = 2. * bcso4oc_col * mac_3im    # effective spherical absorptivity 
r_3im  = 2. * bcso4oc_col * mscb_3im
t_3im  = 1. - a_3im - r_3im

a_2cs  = 2. * bcso4oc_col * mac_2cs    # effective spherical absorptivity 
r_2cs  = 2. * bcso4oc_col * mscb_2cs
t_2cs  = 1. - a_2cs - r_2cs

a_3cs  = 2. * bcso4oc_col * mac_3cs    # effective spherical absorptivity 
r_3cs  = 2. * bcso4oc_col * mscb_3cs
t_3cs  = 1. - a_3cs - r_3cs

alpha = np.zeros(rs.shape)
avgang = np.array([45,40,40,35,45,45,40,45,45,40,40,50,40,40,65])
for i in np.arange(rs.shape[0]):
    alpha[i,:,:] = 1./np.cos(avgang[i]*1./180*np.pi)


#for i in np.arange(a_cs.shape[0]):
#    print (a_cs[i,-1,:] - a_em[i,-1,:])/(r_cs[i,-1,:] - r_em[i,-1,:])

# ignore MSC
#drf_sfc_smod_mon = 2.*f0 * ta * (1. - rs) * 2. * bcso4oc_col * dmac
#drf_toa_smod_mon = 2.*f0 * np.power(ta,2) * 2. * rs * 2. * bcso4oc_col * dmac

# simplified expression
#rf_sfc_em_smod_mon = -1. * f0 * ta * ((1. - rs)*a_em + \
#        r_em*np.power(1-rs,2))
#rf_sfc_cs_smod_mon = -1. * f0 * ta * ((1. - rs)*a_cs + \
#        r_cs*np.power(1-rs,2))
#rf_toa_em_smod_mon = f0 * np.power(ta,2) * (2.*rs*a_em - \
#        r_em*np.power(1-rs,2))
#rf_toa_cs_smod_mon = f0 * np.power(ta,2) * (2.*rs*a_cs - \
#        r_cs*np.power(1-rs,2))

# modified simplified expression with cos(a)
#rf_sfc_em_smod_mon = -1. * f0 * ta * ((1. - rs)*a_em*alpha + \
#        r_em*(1-rs)*(alpha-rs))
#rf_sfc_cs_smod_mon = -1. * f0 * ta * ((1. - rs)*a_cs*alpha + \
#        r_cs*(1-rs)*(alpha-rs))
#rf_toa_em_smod_mon = f0 * np.power(ta,2) * (rs*a_em*(1+alpha) - \
#        r_em*np.power(rs,2) - r_em*alpha + r_em*(1+alpha)*rs)
#rf_toa_cs_smod_mon = f0 * np.power(ta,2) * (rs*a_cs*(1+alpha) - \
#        r_cs*np.power(rs,2) - r_cs*alpha + r_cs*(1+alpha)*rs)
#
#
#drf_sfc_smod_mon = rf_sfc_cs_smod_mon - rf_sfc_em_smod_mon
#drf_toa_smod_mon = rf_toa_cs_smod_mon - rf_toa_em_smod_mon

# full expression
#rf_sfc_em_smod_mon = f0 * ta * ((1. - rs) * (t_em/(1. - rs*r_em) - 1.))
#rf_sfc_cs_smod_mon = f0 * ta * ((1. - rs) * (t_cs/(1. - rs*r_cs) - 1.))
#rf_toa_em_smod_mon = f0 * np.power(ta,2) * (rs - r_em - np.power(t_em,2)*rs/(1. - rs*r_em))
#rf_toa_cs_smod_mon = f0 * np.power(ta,2) * (rs - r_cs - np.power(t_cs,2)*rs/(1. - rs*r_cs))
#drf_sfc_smod_mon = rf_sfc_cs_smod_mon - rf_sfc_em_smod_mon
#drf_toa_smod_mon = rf_toa_cs_smod_mon - rf_toa_em_smod_mon

# modified full expression
t_3em_a  = 1. - a_3em*alpha - r_3em*alpha
t_2im_a  = 1. - a_2im*alpha - r_2im*alpha
t_3im_a  = 1. - a_3im*alpha - r_3im*alpha
t_2cs_a  = 1. - a_2cs*alpha - r_2cs*alpha
t_3cs_a  = 1. - a_3cs*alpha - r_3cs*alpha

rf_sfc_3em_smod_mon = f0 * ta * ((1. - rs) * (t_3em_a/(1. - rs*r_3em) - 1.))
rf_toa_3em_smod_mon = f0 * np.power(ta,2) * (rs - r_3em*alpha - \
        t_3em*t_3em_a*rs/(1. - rs*r_3em))

rf_sfc_2im_smod_mon = f0 * ta * ((1. - rs) * (t_2im_a/(1. - rs*r_2im) - 1.))
rf_toa_2im_smod_mon = f0 * np.power(ta,2) * (rs - r_2im*alpha - \
        t_2im*t_2im_a*rs/(1. - rs*r_2im))

rf_sfc_3im_smod_mon = f0 * ta * ((1. - rs) * (t_3im_a/(1. - rs*r_3im) - 1.))
rf_toa_3im_smod_mon = f0 * np.power(ta,2) * (rs - r_3im*alpha - \
        t_3im*t_3im_a*rs/(1. - rs*r_3im))

rf_sfc_2cs_smod_mon = f0 * ta * ((1. - rs) * (t_2cs_a/(1. - rs*r_2cs) - 1.))
rf_toa_2cs_smod_mon = f0 * np.power(ta,2) * (rs - r_2cs*alpha - \
        t_2cs*t_2cs_a*rs/(1. - rs*r_2cs))

rf_sfc_3cs_smod_mon = f0 * ta * ((1. - rs) * (t_3cs_a/(1. - rs*r_3cs) - 1.))
rf_toa_3cs_smod_mon = f0 * np.power(ta,2) * (rs - r_3cs*alpha - \
        t_3cs*t_3cs_a*rs/(1. - rs*r_3cs))

drf_sfc_2im_smod_mon = rf_sfc_2im_smod_mon - rf_sfc_3em_smod_mon
drf_toa_2im_smod_mon = rf_toa_2im_smod_mon - rf_toa_3em_smod_mon

drf_sfc_3im_smod_mon = rf_sfc_3im_smod_mon - rf_sfc_3em_smod_mon
drf_toa_3im_smod_mon = rf_toa_3im_smod_mon - rf_toa_3em_smod_mon

drf_sfc_2cs_smod_mon = rf_sfc_2cs_smod_mon - rf_sfc_3em_smod_mon
drf_toa_2cs_smod_mon = rf_toa_2cs_smod_mon - rf_toa_3em_smod_mon

drf_sfc_3cs_smod_mon = rf_sfc_3cs_smod_mon - rf_sfc_3em_smod_mon
drf_toa_3cs_smod_mon = rf_toa_3cs_smod_mon - rf_toa_3em_smod_mon


# Average over all months for annual mean
drf_sfc_2im_smod = np.average(drf_sfc_2im_smod_mon,axis=2)
drf_toa_2im_smod = np.average(drf_toa_2im_smod_mon,axis=2)

drf_sfc_3im_smod = np.average(drf_sfc_3im_smod_mon,axis=2)
drf_toa_3im_smod = np.average(drf_toa_3im_smod_mon,axis=2)

drf_sfc_2cs_smod = np.average(drf_sfc_2cs_smod_mon,axis=2)
drf_toa_2cs_smod = np.average(drf_toa_2cs_smod_mon,axis=2)

drf_sfc_3cs_smod = np.average(drf_sfc_3cs_smod_mon,axis=2)
drf_toa_3cs_smod = np.average(drf_toa_3cs_smod_mon,axis=2)

### Two Layer Simplified Mode: Calculate dRF_TOA and dRF_SFC
### assuming equal column density of aerosols

# Define ratio of aerosol column density split between layer 1 (bottom) and
# layer 2 (top)
#acol_r_l1 = 0.5
#acol_r_l2 = 0.5
## calculate t r a for layer 1
#a_em_l1 = acol_r_l1 * a_em
#a_cs_l1 = acol_r_l1 * a_cs
#r_em_l1 = acol_r_l1 * r_em
#r_cs_l1 = acol_r_l1 * r_cs
#t_em_l1 = 1. - a_em_l1 - r_em_l1
#t_cs_l1 = 1. - a_cs_l1 - r_cs_l1
## effective albedo of layer 1 plus surface
#rs_em_l1 = r_em_l1 + np.power(t_em_l1, 2) * rs / (1. - rs * r_em_l1)
#rs_cs_l1 = r_cs_l1 + np.power(t_cs_l1, 2) * rs / (1. - rs * r_cs_l1)
## calculate t r a for layer 2
#a_em_l2 = acol_r_l2 * a_em
#a_cs_l2 = acol_r_l2 * a_cs
#r_em_l2 = acol_r_l2 * r_em
#r_cs_l2 = acol_r_l2 * r_cs
#t_em_l2 = 1. - a_em_l2 - r_em_l2
#t_cs_l2 = 1. - a_cs_l2 - r_cs_l2
#
## calculate TOA forcing
#rf_toa_em_smod_mon2 = f0 * np.power(ta,2) * (rs - r_em_l2 - \
#        np.power(t_em_l2,2)*rs_em_l1/(1. - rs_em_l1*r_em_l2))
#rf_toa_cs_smod_mon2 = f0 * np.power(ta,2) * (rs - r_cs_l2 - \
#        np.power(t_cs_l2,2)*rs_cs_l1/(1. - rs_cs_l1*r_cs_l2))
#drf_toa_smod_mon2 = rf_toa_cs_smod_mon2 - rf_toa_em_smod_mon2
#
## Average over all months for annual mean
#drf_toa_smod2 = np.average(drf_toa_smod_mon2,axis=2)



### Calculate slope from historical data
for region in regions:
    
    i = regions.index(region)
    
    print region
    print 'TOA:'
    print 'EM: '
    print 'Standalone Model:  ',
    print toa_bcsuloc_em[i,:]
    print 'One-Layer SModel:  ',
    print np.average(rf_toa_3em_smod_mon[i,:,:],axis=1)
    print 'Diff:  ',
    print (toa_bcsuloc_em[i,:]-np.average(rf_toa_3em_smod_mon[i,:,:],axis=1))#/np.average(f0[i,:,:]*tau_3em[i,:,:],axis=1)
    print '2IM-EM:  '
    print drf_toa_2im_smod[i,:]
    print '3IM-EM:  '
    print drf_toa_3im_smod[i,:]
    print '2CS-EM:  '
    print drf_toa_2cs_smod[i,:]
    print '3CS-EM:  '
    print drf_toa_3cs_smod[i,:]


    print 'SFC:'
    print 'EM: '
    print 'Standalone Model:  ',
    print sfc_bcsuloc_em[i,:]
    print 'One-Layer SModel:  ',
    print np.average(rf_sfc_3em_smod_mon[i,:,:],axis=1)
    print 'Diff:  ',
    print (sfc_bcsuloc_em[i,:]-np.average(rf_sfc_3em_smod_mon[i,:,:],axis=1))/np.average(f0[i,:,:]*tau_3em[i,:,:],axis=1)
    print '2IM-EM:  '
    print drf_sfc_2im_smod[i,:]
    print '3IM-EM:  '
    print drf_sfc_3im_smod[i,:]
    print '2CS-EM:  '
    print drf_sfc_2cs_smod[i,:]
    print '3CS-EM:  '
    print drf_sfc_3cs_smod[i,:]

    print '\n\n'

### Write results to ncfile
outfile = Dataset('smodforcing_3spec_omr'+str(omr)+'.nc','w',format='NETCDF4')


outfile.createDimension('region',drf_sfc_2im_smod.shape[0])
outfile.createDimension('year',drf_sfc_2im_smod.shape[1])

region_dim = outfile.createVariable('region',np.int32,('region',))
year_dim   = outfile.createVariable('year',np.int32,('year',))

description = ''
for region, rgid in zip(regions, rgids):
    description += region + ': ' + str(rgid) + ',  '
region_dim.description =description

rf_sfc_3em_crtm = outfile.createVariable('rf_sfc_3em_crtm',np.float32,('region','year'))
rf_toa_3em_crtm = outfile.createVariable('rf_toa_3em_crtm',np.float32,('region','year'))

rf_sfc_3em_smod = outfile.createVariable('rf_sfc_3em',np.float32,('region','year'))
rf_toa_3em_smod = outfile.createVariable('rf_toa_3em',np.float32,('region','year'))
rf_sfc_2im_smod = outfile.createVariable('rf_sfc_2im',np.float32,('region','year'))
rf_toa_2im_smod = outfile.createVariable('rf_toa_2im',np.float32,('region','year'))
rf_sfc_3im_smod = outfile.createVariable('rf_sfc_3im',np.float32,('region','year'))
rf_toa_3im_smod = outfile.createVariable('rf_toa_3im',np.float32,('region','year'))
rf_sfc_2cs_smod = outfile.createVariable('rf_sfc_2cs',np.float32,('region','year'))
rf_toa_2cs_smod = outfile.createVariable('rf_toa_2cs',np.float32,('region','year'))
rf_sfc_3cs_smod = outfile.createVariable('rf_sfc_3cs',np.float32,('region','year'))
rf_toa_3cs_smod = outfile.createVariable('rf_toa_3cs',np.float32,('region','year'))

region_dim[:] = rgids
year_dim[:]   = years[:-1]

rf_toa_3em_crtm[:] = toa_bcsuloc_em
rf_sfc_3em_crtm[:] = sfc_bcsuloc_em
rf_toa_3em_smod[:] = np.average(rf_toa_3em_smod_mon,axis=2)
rf_sfc_3em_smod[:] = np.average(rf_sfc_3em_smod_mon,axis=2)
rf_toa_2im_smod[:] = np.average(rf_toa_2im_smod_mon,axis=2)
rf_sfc_2im_smod[:] = np.average(rf_sfc_2im_smod_mon,axis=2)
rf_toa_3im_smod[:] = np.average(rf_toa_3im_smod_mon,axis=2)
rf_sfc_3im_smod[:] = np.average(rf_sfc_3im_smod_mon,axis=2)
rf_toa_2cs_smod[:] = np.average(rf_toa_2cs_smod_mon,axis=2)
rf_sfc_2cs_smod[:] = np.average(rf_sfc_2cs_smod_mon,axis=2)
rf_toa_3cs_smod[:] = np.average(rf_toa_3cs_smod_mon,axis=2)
rf_sfc_3cs_smod[:] = np.average(rf_sfc_3cs_smod_mon,axis=2)

outfile.close()


