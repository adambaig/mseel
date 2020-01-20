# -*- coding: utf-8 -*-
"""
Created on Sat Jan 18 10:41:46 2020

@author: Adam
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import ConvexHull
from sklearn.neighbors import NearestNeighbors

from read_inputs import read_rdvs

F2M = 0.3048
stage = '3H_S20'
grid_spacing= 30
events = read_rdvs(perfs=False)

stages = np.unique([k[:6] for k in events.keys()])

for stage in stages:
    print(stage)
    stage_events =  {k:v for k,v in events.items() if stage in k}
    eastings = [v['easting'] for v in stage_events.values()]
    northings = [v['northing'] for v in stage_events.values()]
    elevations = [v['elevation'] for v in stage_events.values()]
    apparent_stress = [v["SP_ENERGY"]/v["SP_POTENCY"] for v in stage_events.values()]
    fig, ax = plt.subplots(1)
    ax.set_aspect('equal')
    ax.plot(eastings, northings,'.')
    plt.show()
    east_grid = np.arange(min(eastings), max(eastings), grid_spacing)
    north_grid = np.arange(min(northings), max(northings), grid_spacing)
    elev_grid = np.arange(min(elevations), max(elevations), grid_spacing)
    
    #for i_east,east in enumerate(east_grid):
    #    for i_north, north in enumerate(north_grid):
    #        for i_elev, elev in enumerate(elev_grid):
    #            
                
    m_east, m_north, m_elev = np.median(eastings), np.median(northings), np.median(elevations)
    event_keys = []
    for event in stage_events.keys():
        event_keys.append(event)
    
    xyz = np.vstack([eastings, northings, elevations]).T
    nbrs = NearestNeighbors(n_neighbors=10, algorithm='ball_tree').fit(xyz)
    distances, indices = nbrs.kneighbors(xyz)
    plt.plot(np.average(distances))
    clusters = {}
    for event_indices in indices:
        convex_hull = ConvexHull(xyz[event_indices,:])
        volume= convex_hull.volume*F2M*F2M*F2M
        clusters[event_indices[0]] = {}
        clusters[event_indices[0]]["strain"] = sum([stage_events[event_keys[i]]['SP_POTENCY'] for i in event_indices])/2/volume
        clusters[event_indices[0]]["stress"] = sum([stage_events[event_keys[i]]['SP_ENERGY'] for i in event_indices])/clusters[event_indices[0]]["strain"]
        clusters[event_indices[0]]["stiffness"] = clusters[event_indices[0]]["stress"]/clusters[event_indices[0]]["strain"]
        clusters[event_indices[0]]["easting"] = np.median([stage_events[event_keys[i]]['easting'] for i in event_indices])
        clusters[event_indices[0]]["northing"] = np.median([stage_events[event_keys[i]]['northing'] for i in event_indices])
        clusters[event_indices[0]]["elevation"] = np.median([stage_events[event_keys[i]]['elevation'] for i in event_indices])
        clusters[event_indices[0]]["easting"] = np.median([stage_events[event_keys[i]]['easting'] for i in event_indices])
    fig, (ax1,ax2,ax3,ax4) = plt.subplots(1,4,figsize=[20,8],sharex=True, sharey=True)
    for ax in ax1,ax2:
        ax.set_aspect('equal')
        ax.set_facecolor('0.97')
    c_easting= [v['easting'] for v in clusters.values()]
    c_northing= [v['northing'] for v in clusters.values()]
    stiffness = [v["stiffness"] for v in clusters.values()]
    stress = [v["stress"] for v in clusters.values()]
    strain = [v["strain"] for v in clusters.values()]
    ax1.scatter(eastings, northings, c=np.log10(apparent_stress), edgecolor='k', cmap="nipy_spectral")
    ax2.scatter(c_easting, c_northing, c=np.log10(stress), edgecolor='k', cmap="magma")
    ax3.scatter(c_easting, c_northing, c=np.log10(strain), edgecolor='k', cmap="viridis")   
    ax4.scatter(c_easting, c_northing, c=np.log10(stiffness), edgecolor='k', cmap="nipy_spectral") 
    fig.savefig('figures//'+stage+'.png')