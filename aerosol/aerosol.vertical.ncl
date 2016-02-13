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
;region = "China"
;ltlon = 100
;rtlon = 120
;shlat = 25
;nhlat = 50
;shlatb= shlat-1.
;nhlatb= nhlat

;region = "WarmPool"
;ltlon = 90
;rtlon = 130
;shlat = -10
;nhlat = 10
;shlatb = shlat - 1.
;nhlatb = nhlat + 1.

;region = "W-Africa"
;ltlon = 0
;rtlon = 30
;shlat = -15
;nhlat = 10
;shlatb = shlat - 1
;nhlatb = nhlat 

region = "Global Mean"
ltlon = 0
rtlon = 360
shlat = -90
nhlat = 90
shlatb = -90
nhlatb = 90

fact         = 1e9           ;kg -> ug
re           = 6.37122e+6    ; m
pi           = 4.*atan(1.)
d2r          = pi/180.
navo         = 6.02214e23    ; molec/mole
grav         = 9.8           ; m/s2
mw_air       = 28.966        ; g/mole
mw_o3        = 48.           ; g/mole
du_fac       = 21.4
kboltz       = 1.38044e-16 ; erg/K -> g.cm2/s2/K
R            = 287           ; J/(K kg)


;***********
;Options
;***********

plot_filenm = "./plot/aerosol.vertical."+region+".1999"
titletype = "Sulfate (red), BC (black)     "+region
titleunit = "~F8~m~F21~"+"g m~S~-3~N"
;***********
;Read data
;***********

;name & paths
data_filenm = "./aerosol.climatology.AM3.nc"

;read in
filetmp = addfile(data_filenm,"r")

lat     = filetmp->lat({shlat:nhlat})
tlon     = filetmp->lon
lon     = filetmp->lon({ltlon:rtlon})
pfull   = filetmp->pfull
phalf   = filetmp->phalf
time    = filetmp->time
latb    = filetmp->latb({shlatb:nhlatb})

nlat = dimsizes(lat)
ntlon = dimsizes(tlon)
nlon = dimsizes(lon)
npre = dimsizes(pfull)
ntime = dimsizes(time)

;*******************************
;Create and Calculate Ensemble
;*******************************

so4_anthro = filetmp->so4_anthro(ntime-12:ntime-1,:,{shlat:nhlat},{ltlon:rtlon})
so4_natural = filetmp->so4_natural(ntime-12:ntime-1,:,{shlat:nhlat},{ltlon:rtlon})
bc_tmp = filetmp->black_carbon(ntime-12:ntime-1,:,{shlat:nhlat},{ltlon:rtlon})
oc_tmp = filetmp->organic_carbon(ntime-12:ntime-1,:,{shlat:nhlat},{ltlon:rtlon})
so4_tmp = so4_anthro + so4_natural

t_tmp = filetmp->temp(ntime-12:ntime-1,:,{shlat:nhlat},{ltlon:rtlon}) ; get temperature field

;****************************
;Calculate time mean
;****************************
so4_mean = dim_avg_n_Wrap(so4_tmp, 0)
bc_mean  = dim_avg_n_Wrap(bc_tmp, 0)
oc_mean  = dim_avg_n_Wrap(oc_tmp, 0)
t_mean   = dim_avg_n_Wrap(t_tmp,0)

;******************************
;Calculate Area
;******************************

weight = sin(latb(1:nlat)*d2r) - sin(latb(0:nlat-1)*d2r)
area   = 2. * pi * re^2 * weight / ntlon ; m^2

so4_mean1 = wgt_areaave(so4_mean, area, 1.0, 0)
bc_mean1  = wgt_areaave(bc_mean, area, 1.0, 0)
oc_mean1  = wgt_areaave(oc_mean, area, 1.0, 0)

print("so4 burden = "+sum(so4_mean1)*4*pi*re^2*1e-9+"Tg")
print("bc burden = "+sum(bc_mean1)*4*pi*re^2*1e-9+"Tg")
print("oc burden = "+sum(oc_mean1)*4*pi*re^2*1e-9+"Tg")

so4_col = wgt_areaave(so4_mean, area, 1.0, 0) * fact
bc_col  = wgt_areaave(bc_mean, area, 1.0, 0) * fact
oc_col  = wgt_areaave(oc_mean, area, 1.0, 0) * fact

t   = wgt_areaave(t_mean, area, 1.0, 0)

;******************************
;Calculate height of layers
;******************************
phalf(0) = 1e-8

dp = log( phalf(1:npre) / phalf(0:npre - 1) )

dz = R / grav * t * dp

so4 = so4_col / dz
bc  = bc_col / dz
oc  = oc_col / dz

data = new((/3,dimsizes(pfull)/),float)

data(0,:) = so4
data(1,:) = bc
data(2,:) = oc

data!1 = "pfull"
data&pfull = pfull

printVarSummary(data)

;*********
;Plot
;*********
 wks   = gsn_open_wks ("ps",plot_filenm)                   ; open workstation

 res                   = True                       ; plot mods desired
 res@tiMainString      = titletype             ; add title
 res@trYReverse        = True                       ; reverse Y-axis
 res@xyLineThicknesses = (/3.0,3.0,3.0/)               ; make 2nd lines thicker
 res@xyLineColors      = (/"red","black","blue"/)          ; change line color
 res@xyDashPatterns    = (/0.,0.,0./)

  res@tiYAxisString  = "Pressure Level (hPa)"
  res@tiXAxisString  = titleunit
 res@trXMinF  = 0.
; res@trXMaxF  = 30.


 plot  = gsn_csm_xy (wks,data(:,{200:1000}),pfull({200:1000}),res) ; create plot




end
