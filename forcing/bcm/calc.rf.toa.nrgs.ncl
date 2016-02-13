;********************************
;Calculate TOA RF.
;********************************

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

;*****
begin
;*****

years = (/1860,1890,1910,1930,1950,1970,1990/)

rf_di = new((/15,dimsizes(years)/),double)
rf_bc = new((/15,dimsizes(years)/),double)
rf_int  = new((/15,dimsizes(years)/),double)
rf_ext  = new((/15,dimsizes(years)/),double)
rf_sul  = new((/15,dimsizes(years)/),double)
rf_nhyg  = new((/15,dimsizes(years)/),double)
rf_bcsuloc_em = new((/15,dimsizes(years)/),double)

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
data_filenm_nobcsuloc = "./yearly/"+year+"_nobcsuloc.nc"



;read in
filetmp_e = addfile(data_filenm_e,"r")
filetmp_i = addfile(data_filenm_i,"r")
filetmp_ref = addfile(data_filenm_ref,"r")
filetmp_nobc = addfile(data_filenm_nobc,"r")
filetmp_nosul = addfile(data_filenm_nosul,"r")
filetmp_nobcsuloc = addfile(data_filenm_nobcsuloc,"r")

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

olr_e = filetmp_e->olr_clr(0,:,:)
swup_toa_e = filetmp_e->swup_toa_clr(0,:,:)
swdn_toa_e = filetmp_e->swdn_toa_clr(0,:,:)
olr_i = filetmp_i->olr_clr(0,:,:)
swup_toa_i = filetmp_i->swup_toa_clr(0,:,:)
swdn_toa_i = filetmp_i->swdn_toa_clr(0,:,:)
olr_ref = filetmp_ref->olr_clr(0,:,:)
swup_toa_ref = filetmp_ref->swup_toa_clr(0,:,:)
swdn_toa_ref = filetmp_ref->swdn_toa_clr(0,:,:)
olr_nobc = filetmp_nobc->olr_clr(0,:,:)
swup_toa_nobc = filetmp_nobc->swup_toa_clr(0,:,:)
swdn_toa_nobc = filetmp_nobc->swdn_toa_clr(0,:,:)
olr_nosul = filetmp_nosul->olr_clr(0,:,:)
swup_toa_nosul = filetmp_nosul->swup_toa_clr(0,:,:)
swdn_toa_nosul = filetmp_nosul->swdn_toa_clr(0,:,:)
olr_nobcsuloc = filetmp_nobcsuloc->olr_clr(0,:,:)
swup_toa_nobcsuloc = filetmp_nobcsuloc->swup_toa_clr(0,:,:)
swdn_toa_nobcsuloc = filetmp_nobcsuloc->swdn_toa_clr(0,:,:)

rf_e = swdn_toa_e - swup_toa_e - olr_e
rf_i = swdn_toa_i - swup_toa_i - olr_i
rf_ref = swdn_toa_ref - swup_toa_ref - olr_ref
rf_nobc = swdn_toa_nobc - swup_toa_nobc - olr_nobc
rf_nosul = swdn_toa_nosul - swup_toa_nosul - olr_nosul
rf_nobcsuloc = swdn_toa_nobcsuloc - swup_toa_nobcsuloc - olr_nobcsuloc

diff_rf_e = rf_e - rf_ref
diff_rf_i = rf_i - rf_ref
diff_rf_di = diff_rf_i - diff_rf_e
diff_rf_bc = rf_e - rf_nobc
diff_rf_sul = rf_e - rf_nosul
diff_rf_bcsuloc_em = rf_e - rf_nobcsuloc


diff_rf_i!0 = "lat"
diff_rf_i!1 = "lon"
diff_rf_i&lat = lat
diff_rf_i&lon = lon
printVarSummary(diff_rf_i)

copy_VarCoords(diff_rf_i,diff_rf_e)
copy_VarCoords(diff_rf_i,diff_rf_bc)
copy_VarCoords(diff_rf_i,diff_rf_sul)
copy_VarCoords(diff_rf_i,diff_rf_di)
copy_VarCoords(diff_rf_i,diff_rf_bcsuloc_em)

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
diff_rf_bcsuloc_em_mean = wgt_areaave(diff_rf_bcsuloc_em, area, 1.0, 0)

rf_di(0,yr) = diff_rf_di_mean
rf_bc(0,yr) = diff_rf_bc_mean
rf_ext(0,yr) = diff_rf_e_mean
rf_int(0,yr) = diff_rf_i_mean
rf_sul(0,yr) = diff_rf_sul_mean
rf_bcsuloc_em(0,yr) = diff_rf_bcsuloc_em_mean

print("Global: BC+SUL EXT TOA RF = "+diff_rf_e_mean)
print("Global: BC+SUL INT TOA RF = "+diff_rf_i_mean)
print("Global: BC+SUL delta INT TOA RF = "+diff_rf_di_mean)
print("Global: BC EXT TOA RF = "+diff_rf_bc_mean)
print("Global: SUL EXT TOA RF = "+diff_rf_sul_mean)
print("Global: BC+SUL+OC EXT TOA RF = "+diff_rf_bcsuloc_em_mean)

;*****************************************
;Calculate forcing over different regions
;*****************************************
area_globe = new((/nlat,nlon/),float)

do i = 0, nlon - 1
    area_globe(:,i) = doubletofloat(area)
end do

regions = (/"Global","Eastern Asia","Southeastern Asia","Southern Asia","Northern America","Central America","South America","Europe","Northern Africa","Western Africa","Eastern Africa","Southern Africa","Middle Africa","Pacific Warm Pool","Arctic"/)

regionids = (/0, 8, 19, 21, 15, 5, 18, 25, 14, 23, 7, 20, 13, 26, 27/) 

;Read in the mask file
maskfn = "world_regions_2008_"+nlat+"x"+nlon+".nc"
filemask = addfile(maskfn,"r")
rgmask = filemask->mask

do rg = 1, dimsizes(regions) - 1

    region = regions(rg)
    rgid = regionids(rg)

    area_rg = 0.
   
    diff_rf_e_total   = 0.
    diff_rf_i_total   = 0.
    diff_rf_di_total  = 0.
    diff_rf_bc_total  = 0.
    diff_rf_sul_total = 0.
    diff_rf_bcsuloc_em_total = 0.

    do ilat = 0, nlat - 1
        do ilon = 0, nlon - 1
            if (region .eq. "Pacific Warm Pool" ) then
                if ( lat(ilat) .ge. -10 .and. lat(ilat) .le. 10 .and. lon(ilon) .ge. 90 .and. lon(ilon) .le. 130 ) then
                    diff_rf_e_total   =  diff_rf_e_total   + diff_rf_e(ilat, ilon) * area_globe(ilat, ilon)
                    diff_rf_i_total   =  diff_rf_i_total   + diff_rf_i(ilat, ilon) * area_globe(ilat, ilon)
                    diff_rf_di_total  =  diff_rf_di_total  + diff_rf_di(ilat, ilon) * area_globe(ilat, ilon)
                    diff_rf_bc_total  =  diff_rf_bc_total  + diff_rf_bc(ilat, ilon) * area_globe(ilat, ilon)
                    diff_rf_sul_total =  diff_rf_sul_total + diff_rf_sul(ilat, ilon) * area_globe(ilat, ilon)
                    diff_rf_bcsuloc_em_total =  diff_rf_bcsuloc_em_total + diff_rf_bcsuloc_em(ilat, ilon) * area_globe(ilat, ilon)
                    area_rg = area_rg + area_globe(ilat, ilon)
                end if
            else if (region .eq. "Arctic" ) then
                if ( lat(ilat) .ge. 67 ) then
                    diff_rf_e_total   =  diff_rf_e_total   + diff_rf_e(ilat, ilon) * area_globe(ilat, ilon)
                    diff_rf_i_total   =  diff_rf_i_total   + diff_rf_i(ilat, ilon) * area_globe(ilat, ilon)
                    diff_rf_di_total  =  diff_rf_di_total  + diff_rf_di(ilat, ilon) * area_globe(ilat, ilon)
                    diff_rf_bc_total  =  diff_rf_bc_total  + diff_rf_bc(ilat, ilon) * area_globe(ilat, ilon)
                    diff_rf_sul_total =  diff_rf_sul_total + diff_rf_sul(ilat, ilon) * area_globe(ilat, ilon)
                    diff_rf_bcsuloc_em_total =  diff_rf_bcsuloc_em_total + diff_rf_bcsuloc_em(ilat, ilon) * area_globe(ilat, ilon)
                    area_rg = area_rg + area_globe(ilat, ilon)
                end if
            else if ( region .ne. "Europe" ) then 
                if ( rgmask(ilat, ilon) .eq. rgid ) then
                    diff_rf_e_total   =  diff_rf_e_total   + diff_rf_e(ilat, ilon) * area_globe(ilat, ilon)
                    diff_rf_i_total   =  diff_rf_i_total   + diff_rf_i(ilat, ilon) * area_globe(ilat, ilon)
                    diff_rf_di_total  =  diff_rf_di_total  + diff_rf_di(ilat, ilon) * area_globe(ilat, ilon)
                    diff_rf_bc_total  =  diff_rf_bc_total  + diff_rf_bc(ilat, ilon) * area_globe(ilat, ilon)
                    diff_rf_sul_total =  diff_rf_sul_total + diff_rf_sul(ilat, ilon) * area_globe(ilat, ilon)
                    diff_rf_bcsuloc_em_total =  diff_rf_bcsuloc_em_total + diff_rf_bcsuloc_em(ilat, ilon) * area_globe(ilat, ilon)
                    area_rg = area_rg + area_globe(ilat, ilon)
                end if
            else if ( region .eq. "Europe" ) then
                if ( rgmask(ilat, ilon) .eq. 25 \
                    .or. rgmask(ilat, ilon) .eq. 22 \
                    .or. rgmask(ilat, ilon) .eq. 16 \
                    .or. rgmask(ilat, ilon) .eq. 9 ) then
                    diff_rf_e_total   =  diff_rf_e_total   + diff_rf_e(ilat, ilon) * area_globe(ilat, ilon)
                    diff_rf_i_total   =  diff_rf_i_total   + diff_rf_i(ilat, ilon) * area_globe(ilat, ilon)
                    diff_rf_di_total  =  diff_rf_di_total  + diff_rf_di(ilat, ilon) * area_globe(ilat, ilon)
                    diff_rf_bc_total  =  diff_rf_bc_total  + diff_rf_bc(ilat, ilon) * area_globe(ilat, ilon)
                    diff_rf_sul_total =  diff_rf_sul_total + diff_rf_sul(ilat, ilon) * area_globe(ilat, ilon)
                    diff_rf_bcsuloc_em_total =  diff_rf_bcsuloc_em_total + diff_rf_bcsuloc_em(ilat, ilon) * area_globe(ilat, ilon)
                    area_rg = area_rg + area_globe(ilat, ilon)
                end if
            end if
            end if
            end if
            end if
        end do
    end do

    diff_rf_e_rg   =  diff_rf_e_total  / area_rg 
    diff_rf_i_rg   =  diff_rf_i_total  / area_rg
    diff_rf_di_rg  =  diff_rf_di_total / area_rg
    diff_rf_bc_rg  =  diff_rf_bc_total / area_rg
    diff_rf_sul_rg =  diff_rf_sul_total/ area_rg
    diff_rf_bcsuloc_em_rg =  diff_rf_bcsuloc_em_total/ area_rg

    rf_di(rg,yr) = diff_rf_di_rg
    rf_bc(rg,yr) = diff_rf_bc_rg
    rf_sul(rg,yr) = diff_rf_bc_rg
    rf_int(rg,yr) = diff_rf_i_rg
    rf_ext(rg,yr) = diff_rf_e_rg
    rf_bcsuloc_em(rg,yr) = diff_rf_bcsuloc_em_rg


    print(region+"  BC+SUL EXT TOA RF = "+diff_rf_e_rg)
    print(region+"  BC+SUL INT TOA RF = "+diff_rf_i_rg)
    print(region+"  BC+SUL delta INT TOA RF = "+diff_rf_di_rg)
    print(region+"  BC EXT TOA RF = "+diff_rf_bc_rg)
    print(region+"  SUL EXT TOA RF = "+diff_rf_sul_rg)
    print(region+"  BC+SUL+OC EXT TOA RF = "+diff_rf_bcsuloc_em_rg)

end do
end do

system("rm -f toa_"+emtype+"_"+imtype+"_15rgs.nc")
outfile = addfile("toa_"+emtype+"_"+imtype+"_15rgs.nc","c")

rf_di!0 = "region"
rf_di!1 = "year"
rf_di&year   = years
rf_di&region = regionids

description = ""
do i = 0, dimsizes(regions) - 1
    description = description + regions(i) + ": "+regionids(i)+ ",  "
end do

outfile@description = description

copyatt(rf_bc,rf_di)
copyatt(rf_sul,rf_di)
copyatt(rf_int,rf_di)
copyatt(rf_ext,rf_di)
copyatt(rf_bcsuloc_em,rf_di)

outfile->rf_di=rf_di
outfile->rf_bc=rf_bc
outfile->rf_sul=rf_sul
outfile->rf_int=rf_int
outfile->rf_ext=rf_ext
outfile->rf_bcsuloc_em=rf_bcsuloc_em

end
