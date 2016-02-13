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

fact         = 1e18           ;kg/m3 -> ug/m3
re           = 6.37122e+6    ; m
pi           = 4.*atan(1.)
d2r          = pi/180.
navo         = 6.02214e23    ; molec/mole
grav         = 9.8           ; m/s2
mw_air       = 28.966        ; g/mole
mw_o3        = 48.           ; g/mole
du_fac       = 21.4
kboltz       = 1.38044e-16 ; erg/K -> g.cm2/s2/K

;***********
;Options
;***********

plot_filenm = "./plot/so4."+pres+"hPa.1999"
titletype = "Sulfate    "+pres+"hPa"
titleunit = ":F8:m:F21:"+"g m:S:-3:N"
;***********
;Read data
;***********

;name & paths
data_filenm = "./aerosol.climatology.AM3.nc"

;read in
filetmp = addfile(data_filenm,"r")

lat     = filetmp->lat
lon     = filetmp->lon
pfull   = filetmp->pfull
time    = filetmp->time
latb    = filetmp->latb

nlat = dimsizes(lat)
nlon = dimsizes(lon)
npre = dimsizes(pfull)
ntime = dimsizes(time)

;*******************************
;Create and Calculate Ensemble
;*******************************

so4_anthro = filetmp->so4_anthro(ntime-12:ntime-1,{pres},:,:)
so4_natural = filetmp->so4_natural(ntime-12:ntime-1,{pres},:,:)
bc_tmp = filetmp->black_carbon(ntime-12:ntime-1,{pres},:,:)
oc_tmp = filetmp->organic_carbon(ntime-12:ntime-1,{pres},:,:)
so4_tmp = so4_anthro + so4_natural

;******************************
;Calculate Area
;******************************

weight = sin(latb(1:nlat)*d2r) - sin(latb(0:nlat-1)*d2r)
area   = 2. * pi * re^2 * weight / nlon ; m^2

so4_mean = dim_avg_n_Wrap(so4_tmp, 0)
bc_mean  = dim_avg_n_Wrap(bc_tmp, 0)
oc_mean  = dim_avg_n_Wrap(oc_tmp, 0)

so4 = so4_mean
bc = bc_mean
oc = oc_mean

do i = 0, nlat - 1
    so4(i,:) = so4_mean(i,:)/area(i) * fact
    bc(i,:) = bc_mean(i,:)/area(i) * fact
    oc(i,:) = oc_mean(i,:)/area(i) * fact
end do

so4!0 = "lat"
so4!1 = "lon"
so4&lat = lat
so4&lon = lon
printVarSummary(so4)

copy_VarCoords(so4,bc)
copy_VarCoords(so4,oc)

;****************
;Create Plot
;****************

  wks = gsn_open_wks("ps",plot_filenm)               ; Open panel.ps
  gsn_define_colormap(wks,"so4_23")             ; Blue-Red colormap

  res                      = True
  res@gsnMaximize          = True
;  res@gsnDraw              = False              ; don't draw
;  res@gsnFrame             = False              ; don't advance frame

;  res@cnLevelSelectionMode = "AutomaticLevels"     ; set manual contour levels
  res@cnLevelSelectionMode = "ExplicitLevels"

;    res@cnLevels = (/-1, -0.5, 0, 0.1, 0.3, 0.5, 1.0, 1.5, 2.0, 5.0/)
    res@cnLevels = (/-1, -0.5, 0, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0/)
    res@cnLevels = (/-1, -0.5, 0, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 30.0/)
;color

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
;  res@mpMinLonF        = 60.                   ; select a subregion
;  res@mpMaxLonF        = 140.
;  res@mpMinLatF        =  10. 
;  res@mpMaxLatF        =  55.

  res@gsnLeftString    = titletype
  res@gsnRightString   = titleunit
  res@tiMainString     = ""

  res@tmYROn           = False     ; Turn off right and top tick marks
  res@tmXTOn           = False     ; Turn off right and top tick marks

  res@mpGeophysicalLineThicknessF= 2.      ; double the thickness of geophysical boundaries
  res@mpNationalLineThicknessF= 2.         ; double the thickness of national boundaries
  delete(res@mpOutlineSpecifiers)          ; about to change size of attribute, must delete
  res@mpOutlineSpecifiers = "AllNational"

  plots = gsn_csm_contour_map(wks,,res)



end
