;*************************************
;Calculate the aerosol column density.
;*************************************

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

;*****
begin
;*****

;*****************
;Define constants.
;*****************

fact         = 1e9           ;kg -> ug
re           = 6.37122e+6    ; m
pi           = 4.*atan(1.)
d2r          = pi/180.
navo         = 6.02214e23    ; molec/mole
grav         = 9.8           ; m/s2
mw_air       = 28.966        ; g/mole
mw_o3        = 48.           ; g/mole
du_fac       = 21.4
kboltz       = 1.38044e-16 ; erg/K -> g.cm2/s2/K
R            = 287           ; J/(K kg)



regions = (/"Global","China","WarmPool","W-Africa","E-Africa","S-Africa"/)
years   = (/1860,1890,1910,1930,1950,1970,1990,1999/)

bc_col  = new((/dimsizes(regions),dimsizes(years)/),double)
so4_col = new((/dimsizes(regions),dimsizes(years)/),double)
oc_col  = new((/dimsizes(regions),dimsizes(years)/),double)


do rg = 0, dimsizes(regions) - 1

    region = regions(rg)

;***********************
; Define regions
;***********************

if (region .eq. "China") then
ltlon = 100
rtlon = 120
shlat = 25
nhlat = 50
shlatb= shlat-1
nhlatb= nhlat
end if

if (region .eq. "WarmPool") then
ltlon = 90
rtlon = 130
shlat = -10
nhlat = 10
shlatb = shlat - 1
nhlatb = nhlat + 1
end if

if (region .eq. "W-Africa")
ltlon = -20
rtlon = 20
shlat = 0
nhlat = 20
shlatb = shlat - 1
nhlatb = nhlat + 1
end if

if (region .eq. "E-Africa")
ltlon = 20
rtlon = 50
shlat = -10
nhlat = 10
shlatb = shlat - 1
nhlatb = nhlat + 1
end if

if (region .eq. "S-Africa")
ltlon = 20
rtlon = 40
shlat = -25
nhlat = -10
shlatb = shlat - 1
nhlatb = nhlat + 1
end if

if (region .eq. "Global")
ltlon = 0
rtlon = 360
shlat = -90
nhlat = 90
shlatb = -90
nhlatb = 90
end if

;***********
;Read data
;***********

;name & paths
data_filenm = "./aerosol.climatology.AM3.nc"

;read in
filetmp = addfile(data_filenm,"r")

lat     = filetmp->lat({shlat:nhlat})
tlon    = filetmp->lon
lon     = filetmp->lon({ltlon:rtlon})
pfull   = filetmp->pfull
phalf   = filetmp->phalf
time    = filetmp->time
latb    = filetmp->latb({shlatb:nhlatb})

nlat = dimsizes(lat)
ntlon = dimsizes(tlon)
nlon = dimsizes(lon)
npre = dimsizes(pfull)
ntime = dimsizes(time)

;*******************************
;Create and Calculate Ensemble
;*******************************

so4_anthro = filetmp->so4_anthro(:,:,{shlat:nhlat},{ltlon:rtlon})
so4_natural = filetmp->so4_natural(:,:,{shlat:nhlat},{ltlon:rtlon})
bc_tmp = filetmp->black_carbon(:,:,{shlat:nhlat},{ltlon:rtlon})
oc_tmp = filetmp->organic_carbon(:,:,{shlat:nhlat},{ltlon:rtlon})
so4_tmp = so4_anthro + so4_natural

;******************************
;Calculate Area
;******************************

weight = sin(latb(1:nlat)*d2r) - sin(latb(0:nlat-1)*d2r)
area   = 2. * pi * re^2 * weight / ntlon ; m^2

;****************************
;Calculate time mean
;****************************
do yr = 0, dimsizes(years) - 1
    so4_mean = dim_avg_n_Wrap(so4_tmp(yr:yr*12+11,:,:,:), 0)
    bc_mean  = dim_avg_n_Wrap(bc_tmp(yr:yr*12+11,:,:,:), 0)
    oc_mean  = dim_avg_n_Wrap(oc_tmp(yr:yr*12+11,:,:,:), 0)
    
    so4_col(rg,yr) = sum(wgt_areaave(so4_mean, area, 1.0, 0))
    bc_col(rg,yr)  = sum(wgt_areaave(bc_mean, area, 1.0, 0))
    oc_col(rg,yr)  = sum(wgt_areaave(oc_mean, area, 1.0, 0))

    delete(so4_mean)
    delete(bc_mean)
    delete(oc_mean)
end do

delete(lat)
delete(lon)
delete(latb)
delete(weight)
delete(area)
delete(so4_anthro)
delete(so4_natural)
delete(bc_tmp)
delete(oc_tmp)
delete(so4_tmp)

end do

system("rm -f aero_col.nc")
outfile = addfile("aero_col.nc","c")

so4_col!0 = "region"
so4_col!1 = "year"
so4_col&year   = years

copyatt(bc_col,so4_col)
copyatt(oc_col,so4_col)


outfile->so4_col=so4_col
outfile->bc_col=bc_col
outfile->oc_col=oc_col


end
