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


rgs = "_15rgs"

pltfn = "./plot/dint-toa-sfc-"+imscheme+rgs

if ( imscheme .eq. "imbcm_ndist") then
    plttitle = "Homogeneous Mixing - External Mixing"
else if ( imscheme .eq. "csbcm_ndist") then
    plttitle = "Core-shell Mixing - External Mixing"
end if
end if

;*******************
;Read GFDL RTM data
;*******************

ft_toa = addfile("toa_embcm_"+imscheme+rgs+".nc","r")
ft_sfc = addfile("sfc_embcm_"+imscheme+rgs+".nc","r")

toa_di = ft_toa->rf_di
sfc_di = ft_sfc->rf_di

dm = dimsizes(toa_di)
year = dm(1)
region = dm(0)
rg = region

toa_di_plt = toa_di
sfc_di_plt = -1.*sfc_di


;********************
;Read stats
;********************

fh = addfile("./linreg.smodratios."+imscheme+".nc","r")

regstats = fh->regstats

a = regstats(:,0) ;slope
b = regstats(:,1) ;intercept
r = regstats(:,2) ;r-value

linex = new((/rg,100/),double)
liney = new((/rg,100/),double)


do i = 0, rg-1
    llimx = min(toa_di_plt(i,:))
    ulimx = max(toa_di_plt(i,:))
    if ( imscheme .eq. "imbcm_ndist" ) then
        linex(i,:) = fspan(llimx-0.2,ulimx+0.2,100)
    else if ( imscheme .eq. "csbcm_ndist" ) then
        linex(i,:) = fspan(llimx-0.15,ulimx+0.15,100)
    end if
    end if
    liney(i,:) = linex(i,:) * a(i) + b(i)
end do
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
 

    if ( imscheme .eq. "imbcm_ndist" ) then

      res@trYMinF = -0.2
      res@trXMinF = -0.1
      res@trYMaxF = 3.5
      res@trXMaxF = 2.

    else if ( imscheme .eq. "csbcm_ndist" ) then

      res@trYMinF = -0.2
      res@trXMinF = -0.1
      res@trYMaxF = 3.5
      res@trXMaxF = 1.5

    end if
    end if
 
    plots(i)  = gsn_csm_xy(wks,linex(i,:),liney(i,:),res)

    mkres                       = True                      ; polyline mods desired
    mkres@gsMarkerIndex         = i+2                       ; color of lines
    mkres@gsMarkerThicknessF    = 1.5
    mkres@gsMarkerColor         = colors(i+1)               ; thickness of lines
    mkres@gsMarkerSizeF         = 0.011

    dums(i) = gsn_add_polymarker(wks,plots(i),toa_di_plt(i,:),sfc_di_plt(i,:),mkres)

    ; Add annotation

    str_slope = sprintf("%4.2f",a(i))
    str_r2    = sprintf("%4.3f",r(i)*r(i))
    str       = "~F8~a~F21~ = "+str_slope+"~C~~F10~r~S~2~N~~F21~= "+str_r2

    txres                   = True
    txres@txPerimOn         = True
    txres@txFontHeightF     = 0.05 
    txres@txBackgroundFillColor     = "White"

    txid = gsn_create_text(wks, str, txres)

    amres                   = True
    amres@amParallelPosF    = 0.5
    amres@amOrthogonalPosF  = 0.5
    amres@amJust            = "BottomRight"

    annoid = gsn_add_annotation(plots(i), txid, amres)

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
  gsn_text_ndc(wks,"-~F8~D~F10~DRE~B~sfc~N~" + " ~F21~(~F10~W/m~S~2~N~)",0.085,0.5,txres3); add the common left y-axis label

  txres3@txAngleF      = 0.                                               ; put back to normal orientation 
  gsn_text_ndc(wks,"~F8~D~F10~DRE~B~toa~N~" + " ~F21~(~F10~W/m~S~2~N~)",0.5,0.18,txres3); add the common bottom x-axis label
 
  gsn_panel(wks, plots, (/3,5/),resP)

  frame(wks)

  delete(liney)

end do ;imscheme

end
