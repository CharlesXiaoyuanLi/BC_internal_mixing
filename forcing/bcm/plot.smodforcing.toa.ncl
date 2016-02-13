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

do im = 0, 1

print(im)

if ( im .eq. 0 ) then
    imscheme = "imbcm_ndist"
else if ( im .eq. 1 ) then
    imscheme = "csbcm_ndist"
end if
end if

;emtag = "on"
;emtag = ""

rgs = "_15rgs"

;if ( emtag .eq. "on" ) then
;    plttitle = "TOA RF: External Mixing of BC and sulfate"
;    pltfn = "./plot/smodforcings_embcm"+rgs
if ( imscheme .eq. "imbcm_ndist") then
    plttitle = "~F8~D~F~DRE~B~toa~N~: Homogeneous Mixing - External Mixing (BC and sulfate)"
    pltfn   = "./plot/smodforcings_"+imscheme+rgs
else if ( imscheme .eq. "csbcm_ndist") then
    plttitle = "~F8~D~F~DRE~B~toa~N~: Core-shell Mixing - External Mixing (BC and sulfate)"
    pltfn   = "./plot/smodforcings_"+imscheme+rgs
else
    print("ERROR: Please enter 'csbcm_ndist' or 'imbcm_ndist'.")
    exit
end if
end if
;end if

;***********
;Read data
;***********

fh = addfile("./smodforcing_"+imscheme+".nc","r")

rf_sfc_ext_smod = fh->rf_sfc_em
rf_toa_ext_smod = fh->rf_toa_em
if( imscheme .eq. "csbcm_ndist" ) then
    rf_sfc_int_smod = fh->rf_sfc_cs
    rf_toa_int_smod = fh->rf_toa_cs
else if( imscheme .eq. "imbcm_ndist" ) then
    rf_sfc_int_smod = fh->rf_sfc_im
    rf_toa_int_smod = fh->rf_toa_im
else
    print("ERROR: Please enter 'csbcm_ndist' or 'imbcm_ndist'.")
    exit
end if 
end if

rf_toa_drf_smod = rf_toa_int_smod - rf_toa_ext_smod

regionids = fh->region
years     = fh->year

rg = dimsizes(regionids)
yr = dimsizes(years)

fh_toa = addfile("./toa_embcm_"+imscheme+rgs+".nc","r")
fh_sfc = addfile("./sfc_embcm_"+imscheme+rgs+".nc","r")

toa_int = fh_toa->rf_int
toa_ext = fh_toa->rf_ext
toa_drf = toa_int - toa_ext

linex = fspan(-5,5,40)
liney = linex * 1.

;****************
;Plot
;****************
regions = (/"Global Mean","Eastern Asia","Southeastern Asia","Southern Asia","Northern America","Central America","South America","Europe","Northern Africa","Western Africa","Eastern Africa","Southern Africa","Middle Africa","Pacific Warm Pool","Arctic"/)

labels = (/"GL","EAs","SEAs","SAs","NAm","CAm","SAm","EU","NAf","WAf","EAf","SAf","MAf","PWP","AR"/)

colors = (/"white","black","firebrick1","darkgreen","coral4", \ 
             "cyan3","darkslateblue","darkgoldenrod","palevioletred4",     \
             "dodgerblue4","darkslategrey","indianred", \
             "chartreuse4","darkorange3", \
             "purple4","darksalmon"/)

  wks  = gsn_open_wks("ps",pltfn)              ; specifies a ps plot
  gsn_define_colormap(wks,colors)

  res                    = True                   ; plot mods desired

  res@gsnDraw            = False
  res@gsnFrame           = False

  res@gsnMaximize        = True

;Plot line first
  res@xyMarkLineModes     = "Lines"

  res@tmXMajorGrid                = True          ; implement x grid 
  res@tmXMajorGridThicknessF      = 1.0           ; 2.0 is default
  res@tmXMajorGridLineDashPattern = 2             ; select short dash lines
        
  res@tmYMajorGrid                = True          ; implement y grid
  res@tmYMajorGridThicknessF      = 1.0           ; 2.0 is default
  res@tmYMajorGridLineDashPattern = 2             ; select short dash lines

  res@tmXBLabelFontHeightF = 0.04
  res@tmYLLabelFontHeightF = 0.04


  plots = new(rg,"graphic")
  dums = new(rg,"graphic")

  do i = 0, rg - 1

    res@gsnLeftString      = regions(i)
    res@gsnLeftStringFontHeightF = 0.05

;    res@tiXAxisString      = "Linear Regression"
;    res@tiXAxisFontHeightF = 0.021
;    res@tiYAxisString      = "Simplified Model"
;    res@tiYAxisFontHeightF = 0.021
;
    res@trYMinF = -0.1
    res@trXMinF = -0.1
    res@trYMaxF = 1.1
    res@trXMaxF = 1.1
 
    plots(i)  = gsn_csm_xy(wks,linex,liney,res)

    mkres                       = True                      ; polyline mods desired
    mkres@gsMarkerIndex         = i+2                       ; color of lines
    mkres@gsMarkerThicknessF    = 1.5
    mkres@gsMarkerColor         = colors(i+1)               ; thickness of lines
    mkres@gsMarkerSizeF         = 0.011

;    if ( emtag .eq. "on" ) then
;        dums(i) = gsn_add_polymarker(wks,plots(i),toa_ext(i,:),rf_toa_ext_smod(i,:),mkres)
;    else
        dums(i) = gsn_add_polymarker(wks,plots(i),toa_drf(i,:),rf_toa_drf_smod(i,:),mkres)
;    end if

  end do


;*******************
; Draw panel plot  *
;*******************
  resP = True

  resP@txString = plttitle
  resP@txFontHeightF = 0.015

  resP@gsnMaximize = True
  resP@gsnFrame    = False

  resP@gsnPanelXWhiteSpacePercent = 2; set a bit of extra white space between panels in the 
  resP@gsnPanelYWhiteSpacePercent = 2; x and y directions

  resP@gsnPanelLeft  = 0.1; shrink panel plot so that there is extra room for the left y-axis label
;  resP@gsnPanelRight = 0.9; shrink panel plot so that there is extra room for the right y-axis label

  txres3               = True             
  txres3@txAngleF      = 90.          ; Rotate string clockwise 
  txres3@txFontHeightF = 0.015
  gsn_text_ndc(wks,"SRTM" + " (W/m~S~2~N~)",0.085,0.5,txres3); add the common left y-axis label

  txres3@txAngleF      = 0.                                               ; put back to normal orientation 
  gsn_text_ndc(wks,"GFDL RTM" + " (W/m~S~2~N~)",0.5,0.18,txres3); add the common bottom x-axis label
 
  gsn_panel(wks, plots, (/3,5/),resP)

  frame(wks)

end do ; imscheme

end
