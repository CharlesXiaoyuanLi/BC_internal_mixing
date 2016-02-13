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
region = "S-Africa"
season = (/"DJF", "MAM", "JJA", "SON"/)

if (region .eq. "China") then
ltlon = 100
rtlon = 120
shlat = 25
nhlat = 50
shlatb= shlat-1.
nhlatb= nhlat
end if

if (region .eq. "WarmPool") then
ltlon = 90
rtlon = 130
shlat = -10
nhlat = 10
shlatb = shlat - 1.
nhlatb = nhlat + 1.
end if

if (region .eq. "W-Africa")
ltlon = 340
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

plot_filenm = "./plot/aerosol.vertical.season."+region+".1999"
titleunit = ":F8:m:F21:"+"g m:S:-3:N"
;***********
;Read data
;***********

;name & paths
data_filenm = "./aerosol.climatology.AM3.nc"

;read in
filetmp = addfile(data_filenm,"r")

lat     = filetmp->lat({shlat:nhlat})
tlon     = filetmp->lon

if (region .eq. "W-Africa") then
    lon1     = filetmp->lon({ltlon:360})
    lon2     = filetmp->lon({0:rtlon})
    lon_tmp  = (/lon1,lon2/)
    lon      = ndtooned(lon_tmp)
    nlonm    = dimsizes(lon1)
else
    lon     = filetmp->lon({ltlon:rtlon})
end if

pfull   = filetmp->pfull
phalf   = filetmp->phalf
time    = filetmp->time
latb    = filetmp->latb({shlatb:nhlatb})

nlat = dimsizes(lat)
ntlon = dimsizes(tlon)
nlon = dimsizes(lon)
npre = dimsizes(pfull)
ntime = dimsizes(time)

print(nlon)

;*******************************
;Create and Calculate Ensemble
;*******************************

if (region .eq. "W-Africa") then
    so4_anthro = new((/12,npre,nlat,nlon/),float)
    so4_natural= new((/12,npre,nlat,nlon/),float)
    bc_tmp     = new((/12,npre,nlat,nlon/),float)
    oc_tmp     = new((/12,npre,nlat,nlon/),float)
    t_tmp      = new((/12,npre,nlat,nlon/),float)

    so4_anthro(:,:,:,0:nlonm-1) = filetmp->so4_anthro(ntime-12:ntime-1,:,{shlat:nhlat},{ltlon:360})
    so4_anthro(:,:,:,nlonm:) = filetmp->so4_anthro(ntime-12:ntime-1,:,{shlat:nhlat},{0:rtlon})
    so4_natural(:,:,:,0:nlonm-1) = filetmp->so4_natural(ntime-12:ntime-1,:,{shlat:nhlat},{ltlon:360})
    so4_natural(:,:,:,nlonm:) = filetmp->so4_natural(ntime-12:ntime-1,:,{shlat:nhlat},{0:rtlon})
    bc_tmp(:,:,:,0:nlonm-1) = filetmp->black_carbon(ntime-12:ntime-1,:,{shlat:nhlat},{ltlon:360})
    bc_tmp(:,:,:,nlonm:) = filetmp->black_carbon(ntime-12:ntime-1,:,{shlat:nhlat},{0:rtlon})
    oc_tmp(:,:,:,0:nlonm-1) = filetmp->organic_carbon(ntime-12:ntime-1,:,{shlat:nhlat},{ltlon:360})
    oc_tmp(:,:,:,nlonm:) = filetmp->organic_carbon(ntime-12:ntime-1,:,{shlat:nhlat},{0:rtlon})
    t_tmp(:,:,:,0:nlonm-1) = filetmp->temp(ntime-12:ntime-1,:,{shlat:nhlat},{ltlon:360})
    t_tmp(:,:,:,nlonm:) = filetmp->temp(ntime-12:ntime-1,:,{shlat:nhlat},{0:rtlon})


else
    so4_anthro = filetmp->so4_anthro(ntime-12:ntime-1,:,{shlat:nhlat},{ltlon:rtlon})
    so4_natural = filetmp->so4_natural(ntime-12:ntime-1,:,{shlat:nhlat},{ltlon:rtlon})
    bc_tmp = filetmp->black_carbon(ntime-12:ntime-1,:,{shlat:nhlat},{ltlon:rtlon})
    oc_tmp = filetmp->organic_carbon(ntime-12:ntime-1,:,{shlat:nhlat},{ltlon:rtlon})
    t_tmp   = filetmp->temp(:,:,{shlat:nhlat},{ltlon:rtlon})
end if

so4_tmp = so4_anthro + so4_natural

;******************************
;Calculate Area
;******************************

weight = sin(latb(1:nlat)*d2r) - sin(latb(0:nlat-1)*d2r)
area   = 2. * pi * re^2 * weight / ntlon ; m^2

so4_mean = new((/4,npre,nlat,nlon/),float)
bc_mean = new((/4,npre,nlat,nlon/),float)
oc_mean = new((/4,npre,nlat,nlon/),float)
t_mean  = new((/4,npre,nlat,nlon/),float)

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
so4_mean(i,:,:,:) = dim_avg_n_Wrap(so4_tmp((/li,mi,ri/),:,:,:), 0)
bc_mean(i,:,:,:)  = dim_avg_n_Wrap(bc_tmp((/li,mi,ri/),:,:,:), 0)
oc_mean(i,:,:,:)  = dim_avg_n_Wrap(oc_tmp((/li,mi,ri/),:,:,:), 0)
t_mean(i,:,:,:)   = dim_avg_n_Wrap(t_tmp((/li,mi,ri/),:,:,:), 0)
end do

so4_am = wgt_areaave(so4_mean, area, 1.0, 0)
bc_am  = wgt_areaave(bc_mean, area, 1.0, 0)
oc_am  = wgt_areaave(oc_mean, area, 1.0, 0)

t   = wgt_areaave(t_mean, area, 1.0, 0)

;******************************
;Calculate height of layers
;******************************
phalf(0) = 1e-8

dp = log( phalf(1:npre) / phalf(0:npre - 1) )

dz = new((/4,npre/),float)

do i = 0, 3
    dz(i,:) = R / grav * t(i,:) * dp
end do

so4 = so4_am / dz * fact
bc  = bc_am / dz * fact
oc  = oc_am / dz * fact


;*************************
;Calculate column density
;*************************

so4_col = new(4,double)
bc_col  = new(4,double)
oc_col  = new(4,double)

phalf(0) = 1e-8

do i = 0, 3
    so4_col(i) = dim_sum(so4_am(i,:)) * 1e6
    bc_col(i)  = dim_sum(bc_am(i,:)) * 1e6
    oc_col(i)  = dim_sum(oc_am(i,:)) * 1e6
end do

print(so4_col+"mg/m2")
print(bc_col+"mg/m2")
print(oc_col+"mg/m2")


;************************
;Combine data
;************************

data = new((/3,4,dimsizes(pfull)/),float)

data(0,:,:) = so4
data(1,:,:) = bc
data(2,:,:) = oc

ratio = new((/2,4,dimsizes(pfull)/),float)

ratio(0,:,:) = so4/bc
ratio(1,:,:) = oc/bc

data!2 = "pfull"
data&pfull = pfull

ratio!2 = "pfull"
ratio&pfull = pfull

;*********
;Plot
;*********
 wks   = gsn_open_wks ("ps",plot_filenm)                   ; open workstation
 plot = new(4,graphic)

 res                   = True                       ; plot mods desired
 res@gsnDraw           = False
 res@gsnFrame          = False
 res@trYReverse        = True                       ; reverse Y-axis
 res@xyLineThicknesses = (/3.0,3.0,3.0/)               ; make 2nd lines thicker
 res@xyLineColors      = (/"red","black","blue"/)          ; change line color
 res@xyDashPatterns    = (/0.,0.,0./)

  res@tiYAxisString  = "Pressure Level (hPa)"
  res@tiXAxisString  = titleunit
  res@trYMinF  = 300
  res@trYMaxF  = 1000
;54,36,100
  res@xyYStyle = "Log"
;tick marks for Y axis
  res@tmYLMode          = "Explicit"              ; explicit labels
  res@tmYLValues        = (/300,400,600,800,1000/)                       ; location of labels 
  res@tmYLLabels        = ""+res@tmYLValues                ; labels themselves
  res@tmYLMinorValues = ispan(300,1000,100)
  ;res@tmLabelAutoStride = True                    ; nice stride on labels
;no tick mark on right
  res@tmYROn  = False

  res2                   = True                       ; plot mods desired
 res2@gsnDraw           = False
 res2@gsnFrame          = False
 res2@trYReverse        = True                       ; reverse Y-axis
 res2@xyLineThicknesses  = (/2.0, 2.0/)               ; make 2nd lines thicker
 res2@xyLineColors       = (/"red","blue"/)          ; change line color
 res2@xyDashPatterns     = (/11., 11./)

if ( region .eq. "S-Africa" ) then
    res@trXMinF   = 0.
    res@trXMaxF   = 4.
    res2@trXMinF  = 0.
    res2@trXMaxF  = 20.
end if

if ( region .eq. "W-Africa" ) then
    res@trXMinF   = 0.
    res@trXMaxF   = 6.5
    res2@trXMinF  = 0.
    res2@trXMaxF  = 22.
end if

if ( region .eq. "E-Africa" ) then
    res@trXMinF   = 0.
    res@trXMaxF   = 9.
    res2@trXMinF  = 0.
    res2@trXMaxF  = 26.
end if

if ( region .eq. "China" ) then
    res@trXMinF   = 0.
    res@trXMaxF   = 14.
    res2@trXMinF  = 0.
    res2@trXMaxF  = 22.
end if

if ( region .eq. "WarmPool" ) then
    res@trXMinF   = 0.
    res@trXMaxF   = 8.
    res2@trXMinF  = 0.
    res2@trXMaxF  = 24.
end if

if ( region .eq. "Globe" ) then
    res@trXMinF   = 0.
    res@trXMaxF   = 2.
    res2@trXMinF  = 0.
    res2@trXMaxF  = 24.
end if


do i = 0, 3
  res@tiMainOn = True
  res@tiMainString = season(i)
  res@tiMainPosition = "Left"
  plot(i) = gsn_csm_x2y (wks,data(:,i,{300:1000}),ratio(:,i,{300:1000}),pfull({300:1000}),res,res2)
end do

  resp  = True
  resp@txString = "Mass Concentration    "+region
  gsn_panel(wks,plot,(/2,2/),resp)

;***********************************************
; legend resources
;***********************************************
if (region .eq. "China") then

 wksl = gsn_open_wks("ps","./plot/legend.aerosol.vertical.season")

 labels = (/"BC","Sulfate","OM","Sulfate/BC","OM/BC"/)                         ; don't want labels

 lgres                    = True

 lgres@vpWidthF           = 0.7                   ; width of legend (NDC)
 lgres@vpHeightF          = 0.4                   ; height of legend (NDC)


 lgres@lgLineColors       = (/"black","red","blue","red","blue"/)
 lgres@lgDashIndexes      = (/0.,0.,0.,11.,11./)    ; dash indexes
 lgres@lgMonoLineThickness = False
 lgres@lgLineThicknesses  = (/3.0,3.0,3.0,2.0,2.0/)
 lgres@lgLineLabelsOn     = False                ; blank labels

;
; Draw the legend, indicating the number of items, a label for each
; item, and the X, Y position of the legend in NDC coordinates.
;
 gsn_legend_ndc(wksl,5,labels,0.1,0.95,lgres)

end if


end
