# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 14:29:02 2019

@author: Adam
"""

import numpy as np
import matplotlib.pyplot as plt

from read_inputs import read_rdvs


def moment_to_radius(moment,stress_drop):
    return np.cbrt(7.*moment/16./stress_drop)

def radius_to_moment(radius,stress_drop):
    return 7.*stress_drop*radius*radius*radius/16.

events = read_rdvs()
frac_events= {k:v for k,v in events.items() if k[6]=='F'}

energies = np.array([v['SP_ENERGY'] for k,v in frac_events.items()])
moments = np.array([v['SP_MOMENT'] for k,v in frac_events.items()])
radii = np.array([v['SP_RADIUS'] for k,v in frac_events.items()])
stress_drops = np.array([v['SP_STRESS_DROP'] for k,v in frac_events.items()])
snr = np.array([v['QI_SNR'] for k,v in frac_events.items()])
f1,a1 = plt.subplots()
isort = np.argsort(logsnr)
a1.scatter(radii[isort]/3.28,moments[isort],marker='.',c=np.log10(snr[isort]),vmin=0,vmax=2.5 , cmap='Blues')
a1.set_xscale('log')
a1.set_yscale('log')
a1.set_xlabel('source radius (m)')
a1.set_ylabel('seismic moment (N$\cdot$m)')
rad1,rad2 = a1.get_xlim()
mom1,mom2 = a1.get_ylim()
for stress_drop in [1e3,1e4,1e5,1e6]:
    a1.plot([rad1,rad2],[radius_to_moment(rad1,stress_drop),radius_to_moment(rad2,stress_drop)])

a1.set_xlim(rad1,rad2)
a1.set_ylim(mom1,mom2)


f2,a2= plt.subplots()
a2.loglog(snr,stress_drops,'.',alpha=0.05)




