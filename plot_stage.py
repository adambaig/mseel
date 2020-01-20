# -*- coding: utf-8 -*-
"""
Created on Sun Oct  6 14:51:42 2019

@author: Adam
"""
import matplotlib.pyplot as plt
import numpy as np
from read_inputs import read_rdvs, read_wells, read_and_rotate_stage
from datetime import datetime, timedelta


stage_events = read_and_rotate_stage('20','3H')

fig,ax = plt.subplots()

easting = np.array([v['easting'] for k,v in stage_events.items()])
northing = np.array([v['northing'] for k,v in stage_events.items()])
h1 = np.array([v['h1'] for k,v in stage_events.items()])
h2 = np.array([v['h2'] for k,v in stage_events.items()])
elevation = np.array([v['elevation'] for k,v in stage_events.items()])
t_zero = np.array([timedelta(0,0,v['t0']*1.e6) + 
                   datetime.strptime(v['time'],'%d:%m:%y:%H:%M:%S') 
                   for k,v in stage_events.items()])

seconds= np.array([(t-min(t_zero)).seconds for t in t_zero])

ax.set_aspect('equal')
ax = plt.scatter(h1,h2,c=seconds,marker='o', edgecolor='k')


