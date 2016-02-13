### This script reads 3-species internal mixing radiative forcing values
### and plot them in boxes.

import numpy as np
from netCDF4 import Dataset

# read data from ncfile

fn = "./smodforcing_3spec.nc"

fh = Dataset(fn, 'r')

rf_toa_3em_crtm = fh.variables['rf_toa_3em_crtm'][:][:]
rf_toa_3em_smod = fh.variables['rf_toa_3em'][:][:]
rf_toa_2im_smod = fh.variables['rf_toa_2im'][:][:]
rf_toa_3im_smod = fh.variables['rf_toa_3im'][:][:]
rf_toa_2cs_smod = fh.variables['rf_toa_2cs'][:][:]
rf_toa_3cs_smod = fh.variables['rf_toa_3cs'][:][:]

rf_sfc_3em_crtm = fh.variables['rf_sfc_3em_crtm'][:][:]
rf_sfc_3em_smod = fh.variables['rf_sfc_3em'][:][:]
rf_sfc_2im_smod = fh.variables['rf_sfc_2im'][:][:]
rf_sfc_3im_smod = fh.variables['rf_sfc_3im'][:][:]
rf_sfc_2cs_smod = fh.variables['rf_sfc_2cs'][:][:]
rf_sfc_3cs_smod = fh.variables['rf_sfc_3cs'][:][:]

regions = fh.variables['region'][:]

#only calculate 1990 here
pbp = np.zeros((len(regions)-1,5)) #no global mean

pbp[:,0] = rf_toa_3em_crtm[1:,-1]
pbp[:,1] = rf_toa_2im_smod[1:,-1] - rf_toa_3em_smod[1:,-1] + pbp[:,0]
pbp[:,2] = rf_toa_3im_smod[1:,-1] - rf_toa_3em_smod[1:,-1] + pbp[:,0]
pbp[:,3] = rf_toa_2cs_smod[1:,-1] - rf_toa_3em_smod[1:,-1] + pbp[:,0]
pbp[:,4] = rf_toa_3cs_smod[1:,-1] - rf_toa_3em_smod[1:,-1] + pbp[:,0]


import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA

df = pd.DataFrame(pbp, \
        columns=['Mixing\nScheme\nI','Mixing\nScheme\nII','Mixing\nScheme\nIII','Mixing\nScheme\nIV','Mixing\nScheme\nV',])

plt.figure()

bxplt = df.plot(kind='box', return_type='dict', positions=[1,2,3,4,5],whis=1.5)

plt.xlim([0,6])
#plt.ylim([4,27])

# draw minor ticks on y-axis, but turn off on x-axis
plt.minorticks_on()
plt.tick_params(axis='x',which='minor',top='off',bottom='off')
plt.grid(True,
        zorder=0,linestyle=':',linewidth=0.6,which='major',color='k',alpha=0.4,axis='y')
plt.grid(True,
        zorder=0,linestyle=':',linewidth=0.3,which='minor',color='k',alpha=0.4,axis='y')

# configure the color of each boxplot
for i in range(5):
    if i == 0:
        bxplt['boxes'][i].set(color='k',linewidth=1.5)
        bxplt['whiskers'][i*2].set(color='k',linewidth=1.5)
        bxplt['whiskers'][i*2+1].set(color='k',linewidth=1.5)
        bxplt['caps'][i*2].set(color='k',linewidth=1.5)
        bxplt['caps'][i*2+1].set(color='k',linewidth=1.5)
        bxplt['medians'][i].set(color='k',linewidth=1.5)

    elif i == 1:
        bxplt['boxes'][i].set(color='orangered',linewidth=1.5)
        bxplt['whiskers'][i*2].set(color='orangered',linewidth=1.5)
        bxplt['whiskers'][i*2+1].set(color='orangered',linewidth=1.5)
        bxplt['caps'][i*2].set(color='orangered',linewidth=1.5)
        bxplt['caps'][i*2+1].set(color='orangered',linewidth=1.5)
        bxplt['medians'][i].set(color='orangered',linewidth=1.5)

    elif i == 2:
        bxplt['boxes'][i].set(color='maroon',linewidth=1.5)
        bxplt['whiskers'][i*2].set(color='maroon',linewidth=1.5)
        bxplt['whiskers'][i*2+1].set(color='maroon',linewidth=1.5)
        bxplt['caps'][i*2].set(color='maroon',linewidth=1.5)
        bxplt['caps'][i*2+1].set(color='maroon',linewidth=1.5)
        bxplt['medians'][i].set(color='maroon',linewidth=1.5) 

    elif i == 3:
        bxplt['boxes'][i].set(color='c',linewidth=1.5)
        bxplt['whiskers'][i*2].set(color='c',linewidth=1.5)
        bxplt['whiskers'][i*2+1].set(color='c',linewidth=1.5)
        bxplt['caps'][i*2].set(color='c',linewidth=1.5)
        bxplt['caps'][i*2+1].set(color='c',linewidth=1.5)
        bxplt['medians'][i].set(color='c',linewidth=1.5)

    else:
        bxplt['boxes'][i].set(color='teal',linewidth=1.5)
        bxplt['whiskers'][i*2].set(color='teal',linewidth=1.5)
        bxplt['whiskers'][i*2+1].set(color='teal',linewidth=1.5)
        bxplt['caps'][i*2].set(color='teal',linewidth=1.5)
        bxplt['caps'][i*2+1].set(color='teal',linewidth=1.5)
        bxplt['medians'][i].set(color='teal',linewidth=1.5)

for flier in bxplt['fliers']:
    i = bxplt['fliers'].index(flier)
    bxplt['fliers'][i].set(mec='maroon',mew=.8)

# show value of the medians
for line in bxplt['medians']:
    x,y = line.get_xydata()[0] # left point of the median line
    plt.text(x+0.25,y+0.05,'%.1f' % y,
            horizontalalignment='center',verticalalignment='bottom',fontsize=9)

plt.title('$RF_{toa}$:  BC + Sulfate + OC')

# add y-axis label
plt.ylabel('$W / m^{2}$')

plt.draw()
plt.savefig('3spec_toa_rf.ps',format='ps')
