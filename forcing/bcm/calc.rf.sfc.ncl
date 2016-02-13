;********************************
;Calculate Surface RF.
;********************************

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

;*****
begin
;*****

years = (/1860,1890,1910,1930,1950,1970,1990/)

rf_di = new((/6,dimsizes(years)/),double)
rf_bc = new((/6,dimsizes(years)/),double)
rf_int  = new((/6,dimsizes(years)/),double)
rf_ext  = new((/6,dimsizes(years)/),double)
rf_sul  = new((/6,dimsizes(years)/),double)

;***********************
; Define some constants
;***********************

fact         = 1e6           ;kg/m3 -> ug/m3
re           = 6.37122e+6    ; m
pi           = 4.*atan(1.)
d2r          = pi/180.
navo         = 6.02214e23    ; molec/mole
grav         = 9.8           ; m/s2
mw_air       = 28.966        ; g/mole
mw_o3        = 48.           ; g/mole
du_fac       = 21.4
kboltz       = 1.38044e-16 ; erg/K -> g.cm2/s2/K

do yr = 0, dimsizes(years) - 1

year = years(yr)

;***********
;Read data
;***********

;name & paths
data_filenm_e     = "./yearly/"+year+"_"+emtype+".nc"
data_filenm_i     = "./yearly/"+year+"_"+imtype+".nc"
data_filenm_ref   = "./yearly/"+year+"_nobcsul.nc"
data_filenm_nobc  = "./yearly/"+year+"_nobc.nc"
data_filenm_nosul = "./yearly/"+year+"_nosul.nc"

;read in
filetmp_e = addfile(data_filenm_e,"r")
filetmp_i = addfile(data_filenm_i,"r")
filetmp_ref = addfile(data_filenm_ref,"r")
filetmp_nobc = addfile(data_filenm_nobc,"r")
filetmp_nosul = addfile(data_filenm_nosul,"r")

lat     = filetmp_e->lat
lon     = filetmp_e->lon
time    = filetmp_e->time
latb    = filetmp_e->latb

nlat = dimsizes(lat)
nlon = dimsizes(lon)
ntime = dimsizes(time)

;*******************************
;Create and Calculate Ensemble
;*******************************

lwup_sfc_e = filetmp_e->lwup_sfc_clr(0,:,:)
lwdn_sfc_e = filetmp_e->lwdn_sfc_clr(0,:,:)
swup_sfc_e = filetmp_e->swup_sfc_clr(0,:,:)
swdn_sfc_e = filetmp_e->swdn_sfc_clr(0,:,:)
lwup_sfc_i = filetmp_i->lwup_sfc_clr(0,:,:)
lwdn_sfc_i = filetmp_i->lwdn_sfc_clr(0,:,:)
swup_sfc_i = filetmp_i->swup_sfc_clr(0,:,:)
swdn_sfc_i = filetmp_i->swdn_sfc_clr(0,:,:)
lwup_sfc_ref = filetmp_ref->lwup_sfc_clr(0,:,:)
lwdn_sfc_ref = filetmp_ref->lwdn_sfc_clr(0,:,:)
swup_sfc_ref = filetmp_ref->swup_sfc_clr(0,:,:)
swdn_sfc_ref = filetmp_ref->swdn_sfc_clr(0,:,:)
lwup_sfc_nobc = filetmp_nobc->lwup_sfc_clr(0,:,:)
lwdn_sfc_nobc = filetmp_nobc->lwdn_sfc_clr(0,:,:)
swup_sfc_nobc = filetmp_nobc->swup_sfc_clr(0,:,:)
swdn_sfc_nobc = filetmp_nobc->swdn_sfc_clr(0,:,:)
lwup_sfc_nosul = filetmp_nosul->lwup_sfc_clr(0,:,:)
lwdn_sfc_nosul = filetmp_nosul->lwdn_sfc_clr(0,:,:)
swup_sfc_nosul = filetmp_nosul->swup_sfc_clr(0,:,:)
swdn_sfc_nosul = filetmp_nosul->swdn_sfc_clr(0,:,:)

rf_e = swdn_sfc_e - swup_sfc_e - lwup_sfc_e + lwdn_sfc_e
rf_i = swdn_sfc_i - swup_sfc_i - lwup_sfc_i + lwdn_sfc_i
rf_ref = swdn_sfc_ref - swup_sfc_ref - lwup_sfc_ref + lwdn_sfc_ref
rf_nobc = swdn_sfc_nobc - swup_sfc_nobc - lwup_sfc_nobc + lwdn_sfc_nobc
rf_nosul = swdn_sfc_nosul - swup_sfc_nosul - lwup_sfc_nosul + lwdn_sfc_nosul

diff_rf_e = rf_e - rf_ref
diff_rf_i = rf_i - rf_ref
diff_rf_di = diff_rf_i - diff_rf_e
diff_rf_bc = rf_e - rf_nobc
diff_rf_sul = rf_e - rf_nosul


diff_rf_i!0 = "lat"
diff_rf_i!1 = "lon"
diff_rf_i&lat = lat
diff_rf_i&lon = lon
printVarSummary(diff_rf_i)

copy_VarCoords(diff_rf_i,diff_rf_e)
copy_VarCoords(diff_rf_i,diff_rf_bc)
copy_VarCoords(diff_rf_i,diff_rf_sul)
copy_VarCoords(diff_rf_i,diff_rf_di)

;**********************
;Calculate Global Mean
;**********************
weight = sin(latb(1:nlat)*d2r) - sin(latb(0:nlat-1)*d2r)
area   = 2. * pi * re^2 * weight / nlon ; m^2

diff_rf_e_mean = wgt_areaave(diff_rf_e, area, 1.0, 0)
diff_rf_i_mean = wgt_areaave(diff_rf_i, area, 1.0, 0)
diff_rf_di_mean = wgt_areaave(diff_rf_di, area, 1.0, 0)
diff_rf_bc_mean = wgt_areaave(diff_rf_bc, area, 1.0, 0)
diff_rf_sul_mean = wgt_areaave(diff_rf_sul, area, 1.0, 0)

rf_di(0,yr) = diff_rf_di_mean
rf_bc(0,yr) = diff_rf_bc_mean
rf_ext(0,yr) = diff_rf_e_mean
rf_int(0,yr) = diff_rf_i_mean
rf_sul(0,yr) = diff_rf_sul_mean

print("ALL EXT Surface RF = "+diff_rf_e_mean)
print("ALL INT Surface RF = "+diff_rf_i_mean)
print("ALL delta INT Surface RF = "+diff_rf_di_mean)
print("BC EXT Surface RF = "+diff_rf_bc_mean)
print("SUL EXT Surface RF = "+diff_rf_sul_mean)

regions = (/"Global","China","WarmPool","W-Africa","E-Africa","S-Africa"/)

do rg = 1, dimsizes(regions) - 1

region = regions(rg)


if (region .eq. "China") then
ltlon = 100
rtlon = 120
shlat = 25
nhlat = 50
shlatb = shlat - 1.
nhlatb = nhlat + 1.
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

latb_cn = latb({shlatb:nhlatb})
nlat_cn = dimsizes(lat({shlat:nhlat}))

weight_cn = sin(latb_cn(1:nlat_cn)*d2r) - sin(latb_cn(0:nlat_cn-1)*d2r)
area_cn   = 2. * pi * re^2 * weight_cn / nlon ; m^2

diff_rf_e_cn = wgt_areaave(diff_rf_e({shlat:nhlat},{ltlon:rtlon}), area_cn, 1.0, 0)
diff_rf_i_cn = wgt_areaave(diff_rf_i({shlat:nhlat},{ltlon:rtlon}), area_cn, 1.0, 0)
diff_rf_di_cn = wgt_areaave(diff_rf_di({shlat:nhlat},{ltlon:rtlon}), area_cn, 1.0, 0)
diff_rf_bc_cn = wgt_areaave(diff_rf_bc({shlat:nhlat},{ltlon:rtlon}), area_cn, 1.0, 0)
diff_rf_sul_cn = wgt_areaave(diff_rf_sul({shlat:nhlat},{ltlon:rtlon}), area_cn, 1.0, 0)

rf_di(rg,yr) = diff_rf_di_cn
rf_bc(rg,yr) = diff_rf_bc_cn
rf_sul(rg,yr) = diff_rf_bc_cn
rf_int(rg,yr) = diff_rf_i_cn
rf_ext(rg,yr) = diff_rf_e_cn


print(region+"  ALL EXT Surface RF = "+diff_rf_e_cn)
print(region+"  ALL INT Surface RF = "+diff_rf_i_cn)
print(region+"  ALL delta INT Surface RF = "+diff_rf_di_cn)
print(region+"  BC EXT Surface RF = "+diff_rf_bc_cn)
print(region+"  SUL EXT Surface RF = "+diff_rf_sul_cn)

delete(latb_cn)
delete(weight_cn)
delete(area_cn)
end do
end do

system("rm -f sfc_"+emtype+"_"+imtype+".nc")
outfile = addfile("sfc_"+emtype+"_"+imtype+".nc","c")

rf_di!0 = "region"
rf_di!1 = "year"
rf_di&year   = years

copyatt(rf_bc,rf_di)
copyatt(rf_sul,rf_di)
copyatt(rf_int,rf_di)
copyatt(rf_ext,rf_di)


outfile->rf_di=rf_di
outfile->rf_bc=rf_bc
outfile->rf_sul=rf_sul
outfile->rf_int=rf_int
outfile->rf_ext=rf_ext
end
