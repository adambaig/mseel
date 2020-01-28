# -*- coding: utf-8 -*-
"""
Created on Sat Jan 18 10:41:46 2020

@author: Adam
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import ConvexHull, Delaunay
from sklearn.neighbors import NearestNeighbors
from scipy.interpolate import Rbf
from read_inputs import read_rdvs, read_wells

F2M = 0.3048
stage = '3H_S20'
grid_spacing= 30
events = read_rdvs(perfs=False)

wells= read_wells()

stage_events =  {k:v for k,v in events.items() if v["SP_MAGNITUDE"]>-2}

eastings = [v['easting'] for v in stage_events.values()]
northings = [v['northing'] for v in stage_events.values()]
elevations = [v['elevation'] for v in stage_events.values()]
apparent_stress = [v["SP_ENERGY"]/v["SP_POTENCY"] for v in stage_events.values()]

east_grid = np.arange(min(eastings), max(eastings), grid_spacing)
north_grid = np.arange(min(northings), max(northings), grid_spacing)
elev_grid = np.arange(min(elevations), max(elevations), grid_spacing)
            
m_east, m_north, m_elev = np.median(eastings), np.median(northings), np.median(elevations)
event_keys = []
for event in stage_events.keys():
    event_keys.append(event)

def in_hull(p, hull):
    """
    Test if points in `p` are in `hull`

    `p` should be a `NxK` coordinates of `N` points in `K` dimensions
    `hull` is either a scipy.spatial.Delaunay object or the `MxK` array of the 
    coordinates of `M` points in `K`dimensions for which Delaunay triangulation
    will be computed
    """
    if not isinstance(hull,Delaunay):
        hull = Delaunay(hull)

    return hull.find_simplex(p)>=0

xyz = np.vstack([eastings, northings, elevations]).T
nbrs = NearestNeighbors(n_neighbors=10, algorithm='ball_tree').fit(xyz)
distances, indices = nbrs.kneighbors(xyz)
clusters = {}
for event_indices in indices:
    convex_hull = ConvexHull(xyz[event_indices,:])
    volume= convex_hull.volume*F2M*F2M*F2M
    clusters[event_indices[0]] = {}
    clusters[event_indices[0]]["strain"] = sum([stage_events[event_keys[i]]['SP_POTENCY'] for i in event_indices])/2/volume
    clusters[event_indices[0]]["stress"] = sum([stage_events[event_keys[i]]['SP_ENERGY'] for i in event_indices])/clusters[event_indices[0]]["strain"]
    clusters[event_indices[0]]["stiffness"] = clusters[event_indices[0]]["stress"]/clusters[event_indices[0]]["strain"]
    clusters[event_indices[0]]["easting"] = np.average([stage_events[event_keys[i]]['easting'] for i in event_indices])
    clusters[event_indices[0]]["northing"] = np.average([stage_events[event_keys[i]]['northing'] for i in event_indices])
    clusters[event_indices[0]]["elevation"] = np.average([stage_events[event_keys[i]]['elevation'] for i in event_indices])


for ax in ax1,ax2,ax3,ax4:
    ax.set_aspect('equal')
    ax.set_facecolor('0.97')
x1,x2 = ax1.get_xlim()
y1,y2 = ax1.get_ylim()
c_easting= [v['easting'] for v in clusters.values()]
c_northing= [v['northing'] for v in clusters.values()]
c_elevation = [v['elevation'] for v in clusters.values()]
stiffness = [v["stiffness"] for v in clusters.values()]
stress = [v["stress"] for v in clusters.values()]
strain = [v["strain"] for v in clusters.values()]

point_convex_hull = ConvexHull(xyz)

points = np.vstack([c_easting, c_northing, c_elevation]).T
stress_int = Rbf(c_easting, c_northing,c_elevation, np.log10(stress), function='inverse', smooth=0.01)
strain_int = Rbf(c_easting, c_northing,c_elevation, np.log10(strain), function='inverse', smooth=0.01)
stiffness_int = Rbf(c_easting, c_northing,c_elevation, np.log10(stiffness), function='inverse', smooth=0.01, epsilon=10)

ref_elevation= -6380
grid_spacing = 30
east_grid = np.arange(min(eastings),max(eastings),grid_spacing)
north_grid = np.arange(min(northings), max(northings), grid_spacing)
elevation_grid = np.arange(min(elevations), max(elevations), grid_spacing)

ne_grid, nn_grid,nz_grid = len(east_grid), len(north_grid),len(elevation_grid)
stress_map= np.zeros([ne_grid, nn_grid])
strain_map = np.zeros([ne_grid, nn_grid])
stiffness_map = np.zeros([ne_grid, nn_grid])
for i_east,e_grid in enumerate(east_grid):
    for i_north, n_grid in enumerate(north_grid):
#        for i_elev, z_grid in enumerate(elevation_grid):
        stress_map[i_east,i_north] = stress_int(e_grid,n_grid,ref_elevation)
        strain_map[i_east,i_north] = strain_int(e_grid,n_grid,ref_elevation)        
        stiffness_map[i_east,i_north] = stiffness_int(e_grid,n_grid,ref_elevation)
    
fig1, ax1 = plt.subplots(1,figsize=[10,8])
fig2, ax2 = plt.subplots(1,figsize=[10,8])
fig3, ax3 = plt.subplots(1,figsize=[10,8])
fig4, ax4 = plt.subplots(1,figsize=[10,8])
for ax in [ax1,ax2,ax3,ax4]:
    ax.set_aspect('equal')

ax1.scatter(eastings, northings, c=np.log10(apparent_stress), edgecolor='k', cmap="nipy_spectral",zorder=22)
ax2.pcolor(east_grid, north_grid, stress_map.T, cmap="magma",zorder=22,alpha=0.7)
ax3.pcolor(east_grid, north_grid, strain_map.T,  cmap="viridis",zorder=22,alpha=0.7)   
ax4.pcolor(east_grid, north_grid, stiffness_map.T, cmap="nipy_spectral",zorder=22,alpha=0.7)
x1,x2 = ax1.get_xlim()
y1,y2 = ax1.get_ylim()
for ax in ax1,ax2,ax3,ax4:
    for well in wells.values():
        ax.plot(well['easting'],well['northing'],color='firebrick',lw=3,zorder=11)
        ax.plot(well['easting'],well['northing'],color='0.2',lw=2,zorder=12)
        for well_stage,perf in {k:v for k,v in well.items() if k.isdigit()}.items():
            stage_eastings = [v for k,v in perf.items() if 'easting' in k]
            stage_northings = [v for k,v in perf.items() if 'northing' in k]
            ax.plot(stage_eastings,stage_northings,lw=7,color='k',zorder=15)                      
    ax.set_xlim(x1,x2)
    ax.set_ylim(y1,y2)
ax1.set_title('events')
ax2.set_title('clusters\nstress')
ax3.set_title('clusters\nstrain')
ax4.set_title('clusters\nstiffness')
