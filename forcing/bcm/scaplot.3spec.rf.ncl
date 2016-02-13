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

pltfn   = "./plot/scaplot.3spec.rf.toa.ps"
pltfn2  = "./plot/scaplot.3spec.rf.sfc.ps"

;***********
;Read data
;***********

ft = addfile("./smodforcing_3spec.nc","r")

rf_toa_3em_crtm = ft->rf_toa_3em_crtm
rf_toa_3em = ft->rf_toa_3em
rf_toa_2im = ft->rf_toa_2im
rf_toa_3im = ft->rf_toa_3im
rf_toa_2cs = ft->rf_toa_2cs
rf_toa_3cs = ft->rf_toa_3cs

dm = dimsizes(rf_toa_3em)
year = dm(1)
region = dm(0)

rf_sfc_3em_crtm = ft->rf_sfc_3em_crtm
rf_sfc_3em = ft->rf_sfc_3em
rf_sfc_2im = ft->rf_sfc_2im
rf_sfc_3im = ft->rf_sfc_3im
rf_sfc_2cs = ft->rf_sfc_2cs
rf_sfc_3cs = ft->rf_sfc_3cs

rf_toa = new((/region,5/),double)
rf_sfc = new((/region,5/),double)

rf_toa(:,0) = rf_toa_3em_crtm(:,year-1)
rf_toa(:,1) = rf_toa_2im(:,year-1) - rf_toa_3em(:,year-1) + rf_toa(:,0)
rf_toa(:,2) = rf_toa_3im(:,year-1) - rf_toa_3em(:,year-1) + rf_toa(:,0)
rf_toa(:,3) = rf_toa_2cs(:,year-1) - rf_toa_3em(:,year-1) + rf_toa(:,0)
rf_toa(:,4) = rf_toa_3cs(:,year-1) - rf_toa_3em(:,year-1) + rf_toa(:,0)

rf_sfc(:,0) = rf_sfc_3em_crtm(:,year-1)
rf_sfc(:,1) = rf_sfc_2im(:,year-1) - rf_sfc_3em(:,year-1) + rf_sfc(:,0)
rf_sfc(:,2) = rf_sfc_3im(:,year-1) - rf_sfc_3em(:,year-1) + rf_sfc(:,0)
rf_sfc(:,3) = rf_sfc_2cs(:,year-1) - rf_sfc_3em(:,year-1) + rf_sfc(:,0)
rf_sfc(:,4) = rf_sfc_3cs(:,year-1) - rf_sfc_3em(:,year-1) + rf_sfc(:,0)


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

  plot = new(2,graphic)

;---Create a scatter plot
  res                    = True                   ; plot mods desired
;  res@gsnDraw            = False
;  res@gsnFrame           = False
  res@gsnMaximize        = True
  res@gsnLeftString      = "DRE~B~toa~N~"
  res@gsnRightString     = "BC + Sulfate + OC"
  res@gsnLeftStringFontHeightF = 0.021
  res@gsnRightStringFontHeightF  = 0.018


  res@tiYAxisString      = "W/m~S~2~N~"
  res@tiYAxisFontHeightF = 0.018
;---Set some marker resources
  MarkLineModes = new(region,string,"Markers")
  ;MarkLineModes(region) = "Lines"
  res@xyMarkLineModes     = MarkLineModes
  res@xyMarkerThicknessF = 2.5
  res@xyMarkerColors     = colors(1:region)   ; It's okay to list more than you need here

; Set the marker indexes. There are 17 predefined ones at:
  res@xyMarkers          = ispan(2,region+1,1)  ; Again, you can list more than you need.

  res@xyMarkerSizeF      = 0.018

;  res@xyDashPattern     = 0.                    ; Change type, thickness, 
;  res@xyLineThicknessF  = 1.                  ; and color of line
;  res@xyLineColor       = "Black"

  res@trXMinF = 0
  res@trXMaxF = 6


  ; set X-axis tick labels explicitly
  res@tmXBMode      = "Explicit"
  res@tmXBValues    = (/0,1,2,3,4,5,6/)
  cr = "~C~"
;  res@tmXBLabels    = (/""," Mixing"+cr+"Scheme"+cr+"      I"," Mixing"+cr+"Scheme"+cr+"      II"," Mixing"+cr+"Scheme"+cr+"     III"," Mixing"+cr+"Scheme"+cr+"     IV"," Mixing"+cr+"Scheme"+cr+"      V",""/)
  res@tmXBLabels    = (/"","External"+cr+" Mixing"," Homo."+cr+" Mixing"+cr+"BC+Sul"+cr+"   only","Homo."+cr+"Mixing","Core-shell"+cr+"   Mixing"+cr+"  BC+Sul"+cr+"     only","Core-shell"+cr+"   Mixing",""/)

  res@tmXTOn        = False

  res@tmXMajorGrid                = False          ; implement x grid 
        
  res@tmYMajorGrid                = True          ; implement y grid
  res@tmYMajorGridThicknessF      = 1.0           ; 2.0 is default
  res@tmYMajorGridLineDashPattern = 2             ; select short dash lines
  res@tmYMinorGrid                = True
  res@tmYMinorGridThicknessF      = 0.2
  res@tmYMinorGridLineDashPattern = 2

  res@tmXBLabelFontHeightF = 0.012
  res@tmYLLabelFontHeightF = 0.014

;  res@gsnPaperOrientation = "landscape" 
;  res@gsnMaximize         = True

  res@vpWidthF  = 0.5
  res@vpHeightF = 0.8

  plot(0)  = gsn_csm_xy (wks,ispan(1,5,1),rf_toa,res)

  wks2  = gsn_open_wks("ps",pltfn2)              ; specifies a ps plot
  gsn_define_colormap(wks2,colors)

  res2 = res

  res2@gsnLeftString      = "DRE~B~sfc~N~" 
  res2@trYMinF            = -21.
  res2@trYReverse          = True

;  res@gsnLeftString      = "RF~B~sfc~N~" 
;  res@trYMinF            = -19.

;  delete(res@tiYAxisString)

  plot(1) = gsn_csm_xy (wks2,ispan(1,5,1),rf_sfc,res2)

;  gsn_panel(wks,plot,(/1,2/),False)

exit
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
