;********************************
;Plot the contour of Species.
;********************************

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

;*****
begin
;*****

;***********************
; Define some constants
;***********************
pres = 900

season = (/"DJF", "MAM", "JJA", "SON"/)

fact         = 1e9           ; kg -> ug
re           = 6.37122e+6    ; m
pi           = 4.*atan(1.)
d2r          = pi/180.
navo         = 6.02214e23    ; molec/mole
grav         = 9.8           ; m/s2
mw_air       = 28.966        ; g/mole
mw_o3        = 48.           ; g/mole
du_fac       = 21.4
kboltz       = 1.38044e-16 ; erg/K -> g.cm2/s2/K
R            = 287           ; J / (kg K)

;***********
;Options
;***********

plot_filenm = "./plot_shine/sul.dist.1999"
titleunit = "~F8~m~F~g m~S~-3~N~"
;***********
;Read data
;***********

;name & paths
data_filenm = "./aerosol.climatology.AM3.nc"

;read in
filetmp = addfile(data_filenm,"r")

lat     = filetmp->lat(:)
tlon     = filetmp->lon
lon     = filetmp->lon(:)

pfull   = filetmp->pfull
phalf   = filetmp->phalf
time    = filetmp->time
latb    = filetmp->latb(:)

nlat = dimsizes(lat)
ntlon = dimsizes(tlon)
nlon = dimsizes(lon)
npre = dimsizes(pfull)
ntime = dimsizes(time)

print(nlon)

;*******************************
;Create and Calculate Ensemble
;*******************************

    so4_anthro = filetmp->so4_anthro(ntime-12:ntime-1,{pres},:,:)
    so4_natural = filetmp->so4_natural(ntime-12:ntime-1,{pres},:,:)
    bc_tmp = filetmp->black_carbon(ntime-12:ntime-1,{pres},:,:)
    oc_tmp = filetmp->organic_carbon(ntime-12:ntime-1,{pres},:,:)
    t_tmp   = filetmp->temp(ntime-12:ntime-1,{pres},:,:)

so4_tmp = so4_anthro + so4_natural

;******************************
;Calculate Area
;******************************

weight = sin(latb(1:nlat)*d2r) - sin(latb(0:nlat-1)*d2r)
area   = 2. * pi * re^2 * weight / ntlon ; m^2

so4_mean = new((/4,nlat,nlon/),float)
bc_mean = new((/4,nlat,nlon/),float)
oc_mean = new((/4,nlat,nlon/),float)
t_mean  = new((/4,nlat,nlon/),float)

do i = 0, 3
if ( i .eq. 0 ) then
li = 0
mi = 1
ri = 11
else if (i .eq. 1) then
li = 2
mi = 3
ri = 4
else if (i .eq. 2) then
li = 5
mi = 6
ri = 7
else if (i .eq. 3) then
li = 8
mi = 9
ri = 10
end if
end if
end if
end if
so4_mean(i,:,:) = dim_avg_n_Wrap(so4_tmp((/li,mi,ri/),:,:), 0)
bc_mean(i,:,:)  = dim_avg_n_Wrap(bc_tmp((/li,mi,ri/),:,:), 0)
oc_mean(i,:,:)  = dim_avg_n_Wrap(oc_tmp((/li,mi,ri/),:,:), 0)
t_mean(i,:,:)   = dim_avg_n_Wrap(t_tmp((/li,mi,ri/),:,:), 0)
end do

;******************************
;Calculate height of layers
;******************************
dp = log( phalf(npre-7) / phalf(npre-8) ) ;pres = 900

dz = new((/4,nlat,nlon/),float)

dz = R / grav * t_mean * dp

so4 = so4_mean / dz * fact
bc  = bc_mean / dz * fact
oc  = oc_mean / dz * fact

so4!0 = "season"
so4!1 = "lat"
so4!2 = "lon"
so4&lat = lat
so4&lon = lon
printVarSummary(so4)

copy_VarCoords(so4,bc)
copy_VarCoords(so4,oc)

;****************
;Create Plot
;****************
do i = 0,2,2

  wks = gsn_open_wks("ps",plot_filenm+"-"+season(i))               ; Open panel.ps
  gsn_define_colormap(wks,"so4_23")             ; Blue-Red colormap

  res                      = True
  res@gsnMaximize          = True
;  res@gsnDraw              = False              ; don't draw
;  res@gsnFrame             = False              ; don't advance frame

;  res@cnLevelSelectionMode = "AutomaticLevels"     ; set manual contour levels
  res@cnLevelSelectionMode = "ExplicitLevels"


if ( i .eq. 0 )
    res@cnLevels = (/-1, -0.5, 0, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0/)  ; sulfate jja
else
    res@cnLevels = (/-1, -0.5, 0, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0/)
end if

;   res@cnLevels = (/-1, -0.5, 0, 0.05, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0/)

;  res@cnMinLevelValF       =  0.                ; set min contour level
;  res@cnMaxLevelValF       =  1.                ; set max contour level
;  res@cnLevelSpacingF      =  0.1               ; set contour spacing

  res@cnFillOn         = True
  res@cnLinesOn        = False
  res@cnLineLabelsOn   = False
  res@cnFillDrawOrder  = "Draw"             ; fill and lines before map
  res@lbLabelBarOn     = True
  res@lbOrientation    = "Vertical"

;  res@gsnSpreadColors   = True                 ; Span full colormap
;  res@gsnSpreadColorEnd = -2                   ; Don't include gray

;  res@mpCenterLonF     = 180.                   ; center plot at 180
  res@mpMinLonF        = -30.                   ; select a subregion
  res@mpMaxLonF        = 50.
  res@mpMinLatF        =  -30.
  res@mpMaxLatF        =  30.

  res@gsnLeftString    = season(i)
  res@gsnRightString   = titleunit

    res@tmYROn           = False     ; Turn off right and top tick marks
  res@tmXTOn           = False     ; Turn off right and top tick marks

  res@mpGeophysicalLineThicknessF= 2.      ; double the thickness of geophysical boundaries
  res@mpNationalLineThicknessF= 2.         ; double the thickness of national boundaries
  delete(res@mpOutlineSpecifiers)          ; about to change size of attribute, must delete
  res@mpOutlineSpecifiers = "AllNational"

  plots = gsn_csm_contour_map(wks,so4(i,:,:),res)

end do

end
