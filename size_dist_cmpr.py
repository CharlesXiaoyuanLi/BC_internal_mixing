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

D_pgs = [0.0118, 0.0355, 0.05]

#D_pgs = [0.0118, 0.0118, 0.0118]

sigma_g = 2.

x = np.logspace(-4,1,num=10000)

colors = ('r','b','y')

for D_pg, color in zip(D_pgs,colors):
    n = lognormal_dist_log10(x, D_pg, sigma_g, N)
    plt.plot(x, np.pi/6.*np.power(x,3)*n, linewidth=2, color=color)
    plt.xscale('log')
    plt.axis('tight')
    def lognormal_dist_vol(rad):
        return np.pi/6.*np.power(rad,3)*N/(np.sqrt(2. * np.pi) * rad * np.log(sigma_g)) * np.exp( -1. *
                            np.power(np.log(rad/D_pg),2) / (2. *
                                np.power(np.log(sigma_g),2)) )
    def lognormal_dist_log10_vol(rad):
        return np.pi/6.*np.power(rad,3)*2.303*N/(np.sqrt(2. * np.pi) * np.log(sigma_g)) * np.exp( -1. *
                            np.power(np.log(rad/D_pg),2) / (2. *
                                np.power(np.log(sigma_g),2)) )



    print 'D_pg = '+str(D_pg)
    volpartial = quadrature(lognormal_dist_vol,0.1,1e0,tol=3e-11,maxiter=1000)[0]
    voltotal = quadrature(lognormal_dist_vol,1e-3,1e0,tol=3e-11,maxiter=1000)[0]
    print 'Volume Ratio of D_p from 0.1 to 1 = '+str(volpartial/voltotal)
plt.savefig('lognormal.ps',format='ps')
