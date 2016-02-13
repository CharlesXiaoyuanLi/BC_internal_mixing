;********************************
;Plot RF.
;********************************

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

;*****
begin
;*****

rgs = "_15rgs"

magtag = ""

pltfn   = "./plot/dint-toa-sfc-"+emtype+"_"+imtype+rgs+magtag

ystring = "-~F8~D~F~RF~B~SFC~N~" 
xstring = "~F8~D~F~RF~B~TOA~N~"

;***********
;Read data
;***********

ft_toa = addfile("toa_"+emtype+"_"+imtype+rgs+".nc","r")
ft_sfc = addfile("sfc_"+emtype+"_"+imtype+rgs+".nc","r")

toa_di = ft_toa->rf_di

dm = dimsizes(toa_di)
year = dm(1)
region = dm(0)

sfc_di = ft_sfc->rf_di

toa_di_plt = new((/region,year/),double)
sfc_di_plt = new((/region,year/),double)

toa_di_plt(0:region-1,:) = toa_di
sfc_di_plt(0:region-1,:) = -1.*sfc_di
;toa_di_plt(region,:) = fspan(0,3,year)
;if (imtype .eq. "csbcm") then
;    sfc_di_plt(region,:) = toa_di_plt(region,:)*4.4/1.6
;else if (imtype .eq. "imbcm") then
;    sfc_di_plt(region,:) = toa_di_plt(region,:)*2.2/1.6
;else if (imtype .eq. "imbcm_ndist") then
;    sfc_di_plt(region,:) = toa_di_plt(region,:)*2.8/1.6
;else if (imtype .eq. "csbcm_ndist") then
;    sfc_di_plt(region,:) = toa_di_plt(region,:)*9.5/3.
;else
;    print("ERROR: no imtype matched")
;    exit
;end if
;end if
;end if
;end if


;do i = 0, 5
;toa_di_plt(i,:) = toa_di_plt(i,:) -toa_di_plt(i,0)
;end do
;****************
;Plot
;****************
labels = (/"Global Mean","Eastern Asia","Southeastern Asia","Southern Asia","Northern America","Central America","South America","Europe","Northern Africa","Western Africa","Eastern Africa","Southern Africa","Middle Africa","Pacific Warm Pool","Arctic"/)

labels = (/"GL","EA","SEA","SA","NA","CA","SA","EU","NA","WA","EA","SA","MA","PWP","AR"/)

colors = (/"white","black","firebrick1","darkgreen","coral4", \ 
             "cyan3","darkslateblue","darkgoldenrod","palevioletred4",     \
             "dodgerblue4","darkslategrey","indianred", \
             "chartreuse4","darkorange3", \
             "purple4","darksalmon"/)

  wks  = gsn_open_wks("ps",pltfn)              ; specifies a ps plot
  gsn_define_colormap(wks,colors)

;---Create a scatter plot
  res                    = True                   ; plot mods desired
  res@gsnMaximize        = True
  res@gsnRightString      = "W m~S~-2~N"
  res@gsnRightStringFontHeightF = 0.017
  if (imtype .eq. "imbcm" .or. imtype .eq. "imbcm_ndist") then 
    res@gsnLeftString      = "Homogeneous Mixing - External Mixing"
  else if (imtype .eq. "csbcm" .or. imtype .eq. "csbcm_ndist") then
    res@gsnLeftString      = "Core-shell Mixing - External Mixing" 
  else
    print("ERROR: no imtype matched")
    exit
  end if
  end if
  res@gsnLeftStringFontHeightF = 0.021


  res@tiXAxisString      = xstring
  res@tiXAxisFontHeightF = 0.021
  res@tiYAxisString      = ystring
  res@tiYAxisFontHeightF = 0.021
;---Set some marker resources
  MarkLineModes = new(region,string,"Markers")
  ;MarkLineModes(region) = "Lines"
  res@xyMarkLineModes     = MarkLineModes
  res@xyMarkerThicknessF = 2.5
  res@xyMarkerColors     = colors(1:region)   ; It's okay to list more than you need here

; Set the marker indexes. There are 17 predefined ones at:
  res@xyMarkers          = ispan(2,region+1,1)  ; Again, you can list more than you need.

  res@xyMarkerSizeF      = 0.012
  if ( magtag .ne. "") then
      res@xyMarkerSizeF  = 0.0135
  end if

;  res@xyDashPattern     = 0.                    ; Change type, thickness, 
;  res@xyLineThicknessF  = 1.                  ; and color of line
;  res@xyLineColor       = "Black"

  res@trYMinF = 0.
  res@trXMinF = -0.1
;  if ( imtype .eq. "imbcm" .or. imtype .eq. "csbcm" ) then
;      res@trYMaxF = 6.
;      res@trXMaxF = 1.6
;  else if ( imtype .eq. "imbcm_ndist" .or. imtype .eq. "csbcm_ndist" ) then
;      res@trYMaxF = 10.
;      res@trXMaxF = 3.
;  end if
;  end if

  res@trYMaxF = 3.
  res@trXMaxF = 1.8


  res@tmXMajorGrid                = True          ; implement x grid 
  res@tmXMajorGridThicknessF      = 1.0           ; 2.0 is default
  res@tmXMajorGridLineDashPattern = 2             ; select short dash lines
        
  res@tmYMajorGrid                = True          ; implement y grid
  res@tmYMajorGridThicknessF      = 1.0           ; 2.0 is default
  res@tmYMajorGridLineDashPattern = 2             ; select short dash lines

  res@tmXBLabelFontHeightF = 0.015
  res@tmYLLabelFontHeightF = 0.015

  if (magtag .ne. "") then
     res@tmXBLabelFontHeightF = 0.018
     res@tmYLLabelFontHeightF = 0.018
  end if
 

  plot  = gsn_csm_xy (wks,toa_di_plt,sfc_di_plt,res)

;*************
;Plot legend
;*************

wksl = gsn_open_wks("ps","./plot/dint-legend")

lgres               = True

lgres@vpWidthF      = 0.3
lgres@vpHeightF     = 0.8

lgres@lgMonoItemType        = False
lgres@lgMarkerSizeF         = 0.02
lgres@lgMarkerColors        = res@xyMarkerColors(::-1)      ;reverse the array to have Global on top
lgres@lgItemTypes           = new(region,string,"Markers")
lgres@lgMarkerIndexes       = res@xyMarkers(::-1)
lgres@lgMarkerThicknessF    = res@xyMarkerThicknessF
lgres@lgLabelFontHeightF    = 0.05

gsn_legend_ndc(wksl,region,labels(::-1),0.1,0.95,lgres)

end
