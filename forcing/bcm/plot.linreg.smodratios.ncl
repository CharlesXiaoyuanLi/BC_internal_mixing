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

if (im .eq. 0) then
    imtype = "imbcm_ndist"
    pltcentertitle = "Homogeneous Mixing - External Mixing"
else if (im .eq. 1) then
    imtype = "csbcm_ndist"
    pltcentertitle = "Core-shell Mixing - External Mixing"
end if
end if

rgs = "_15rgs"

pltfn   = "./plot/linreg_smodratios_"+imtype+"_"+rgs

;***********
;Read data
;***********

fh = addfile("./linreg.smodratios."+imtype+".nc","r")

regstats_in   = fh->regstats
smodratios_in = fh->smodratios

regionids = fh->region
years     = fh->year

rg = dimsizes(regionids)
yr = dimsizes(years)

; Create variables for plotting, the extra dimension is for y=x line
regslopes  = new(rg,float)
smodratios = new(rg,float)

smodratios = dim_avg_n(smodratios_in,1)
regslopes = regstats_in(:,0)

linex = fspan(0,10,20)
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
  res@gsnRightString     = "W m~S~-2~N"
  res@gsnRightStringFontHeightF = 0.017
  res@tiMainString      = "|~F8~D~F~DRE~B~SFC~N~| / |~F8~D~F~DRE~B~TOA~N~|"
  res@tiMainFontHeightF = 0.021
  res@gsnLeftString    = pltcentertitle
  res@gsnLeftStringFontHeightF = 0.017


  res@tiXAxisString      = "GFDL RTM"
  res@tiXAxisFontHeightF = 0.021
  res@tiYAxisString      = "SRTM"
  res@tiYAxisFontHeightF = 0.021

;Plot line first
  res@xyMarkLineModes     = "Lines"

  res@trYMinF = 0.
  res@trXMinF = 0.
  res@trYMaxF = 9.5
  res@trXMaxF = 9.5

  res@tmXMajorGrid                = True          ; implement x grid 
  res@tmXMajorGridThicknessF      = 1.0           ; 2.0 is default
  res@tmXMajorGridLineDashPattern = 2             ; select short dash lines
        
  res@tmYMajorGrid                = True          ; implement y grid
  res@tmYMajorGridThicknessF      = 1.0           ; 2.0 is default
  res@tmYMajorGridLineDashPattern = 2             ; select short dash lines

  res@tmXBLabelFontHeightF = 0.015
  res@tmYLLabelFontHeightF = 0.015

  plot  = gsn_csm_xy(wks,linex,liney,res)

; Plot markers

  dums = new(rg,graphic)

  do i = 0, rg - 1

    mkres                       = True                      ; polyline mods desired
    mkres@gsMarkerIndex         = i+2                       ; color of lines
    mkres@gsMarkerThicknessF    = 3.0
    mkres@gsMarkerColor         = colors(i+1)               ; thickness of lines
    mkres@gsMarkerSizeF         = 0.025

    dums(i) = gsn_add_polymarker(wks,plot,regslopes(i),smodratios(i),mkres)
  end do

  draw(plot)
  frame(wks)


;*************
;Plot legend
;*************

wksl = gsn_open_wks("ps","./plot/linreg-legend")

lgres               = True

lgres@vpWidthF      = 0.15
lgres@vpHeightF     = 0.8

lgres@lgMonoItemType        = False
lgres@lgMarkerSizeF         = 0.02
lgres@lgMarkerColors        = colors(rg:1)      ;reverse the array to have Global on top
lgres@lgItemTypes           = new(rg,string,"Markers")
lgres@lgMarkerIndexes       = ispan(rg+1,2,1)
lgres@lgMarkerThicknessF    = 3.0
lgres@lgLabelFontHeightF    = 0.1

gsn_legend_ndc(wksl,rg,labels(::-1),0.1,0.95,lgres)

end do ;imtype

end
