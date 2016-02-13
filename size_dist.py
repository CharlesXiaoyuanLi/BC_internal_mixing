#!/bin/python

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate.quadrature as quadrature

def lognormal_dist(rad, D_pg, sigma_g, N):
    return N/(np.sqrt(2. * np.pi) * rad * np.log(sigma_g)) * np.exp( -1. *
            np.power(np.log(rad/D_pg),2) / (2. * np.power(np.log(sigma_g),2)) )
def lognormal_dist_log10(rad, D_pg, sigma_g, N):
    return 2.303*N/(np.sqrt(2. * np.pi) * np.log(sigma_g)) * np.exp( -1. *
            np.power(np.log(rad/D_pg),2) / (2. * np.power(np.log(sigma_g),2)) )



N = 100

D_pg, sigma_g = 0.0118, 2.

x = np.logspace(-4,1,num=10000)

fig, ax1 = plt.subplots()

n = lognormal_dist(x, D_pg, sigma_g, N)
#lns1 = ax1.plot(x, n, linewidth=1.5, color='b',label='$n(D_p)$')
n10 = lognormal_dist_log10(x, D_pg, sigma_g, N)
lns1 = ax1.plot(x, n10, linewidth=1.5, color='b',label='$n(log_{10}D_p)$')
ax1.set_xscale('log')
ax1.set_xlabel('$\mu m$')
ax1.set_ylabel('$\mu m^{-1}$')
#ax1.axis('tight')

def lognormal_dist_vol(rad):
    return np.pi/6.*np.power(rad,3)*N/(np.sqrt(2. * np.pi) * rad * np.log(sigma_g)) * np.exp( -1. *
                        np.power(np.log(rad/D_pg),2) / (2. *
                            np.power(np.log(sigma_g),2)) )
def lognormal_dist_log10_vol(rad):
    return np.pi/6.*np.power(rad,3)*2.303*N/(np.sqrt(2. * np.pi) * np.log(sigma_g)) * np.exp( -1. *
                        np.power(np.log(rad/D_pg),2) / (2. *
                            np.power(np.log(sigma_g),2)) )

print 'D_pg = '+str(D_pg)
volpartial = quadrature(lognormal_dist_vol,0.2,1e0,tol=3e-11,maxiter=1000)[0]
voltotal = quadrature(lognormal_dist_vol,1e-3,1e0,tol=3e-11,maxiter=1000)[0]
print 'Volume Ratio of D_p from 0.1 to 1 = '+str(volpartial/voltotal)

nv = np.pi/6.*np.power(x,3)*n
nv10 = np.pi/6.*np.power(x,3)*n10
ax2 = ax1.twinx()
#lns2 = ax2.plot(x, nv, linewidth=1.5, color='r',label='$n_v(logD_p)$')
lns2 = ax2.plot(x, nv10, linewidth=1.5, color='r',label='$n_v(log_{10}D_p)$')
xfill = np.logspace(-1,1,num=1000)
#yfill = lognormal_dist_vol(xfill)
yfill = lognormal_dist_log10_vol(xfill)
ax2.fill_between(x=xfill,y1=yfill,y2=0,color='r',facecolor='red',alpha=0.5,linewidth=0)

ax2.set_ylabel('$\mu m^3 \cdot \mu m^{-1}$')

lns = lns1 + lns2
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs)

plt.title('$D_{pg} = 0.0118 \mu m$, $\sigma_g = 2.0 \mu m$',
        loc='left',fontsize=13)
plt.tight_layout()
#plt.savefig('size_dist_log10.pdf',format='pdf')
plt.show()
