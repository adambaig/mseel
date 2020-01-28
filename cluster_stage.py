# -*- coding: utf-8 -*-
"""
Created on Sat Jan 18 10:41:46 2020

@author: Adam
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import ConvexHull
from sklearn.neighbors import NearestNeighbors

from read_inputs import read_rdvs, read_wells

F2M = 0.3048
stage = '3H_S20'
grid_spacing= 30
events = read_rdvs(perfs=False)

wells= read_wells()
stages = np.unique([k[:6] for k in events.keys()])

for stage in stages:
    stage_events =  {k:v for k,v in events.items() if stage in k if v["SP_MAGNITUDE"]>-2}
    if len(stage_events) > 50:
        eastings = [v['easting'] for v in stage_events.values()]
        northings = [v['northing'] for v in stage_events.values()]
        elevations = [v['elevation'] for v in stage_events.values()]
        apparent_stress = [v["SP_ENERGY"]/v["SP_POTENCY"] for v in stage_events.values()]
    
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
        for ax in ax1,ax2,ax3,ax4:
            ax.set_aspect('equal')
            ax.set_facecolor('0.97')
        x1,x2 = ax1.get_xlim()
        y1,y2 = ax1.get_ylim()
        c_easting= [v['easting'] for v in clusters.values()]
        c_northing= [v['northing'] for v in clusters.values()]
        stiffness = [v["stiffness"] for v in clusters.values()]
        stress = [v["stress"] for v in clusters.values()]
        strain = [v["strain"] for v in clusters.values()]
        ax1.scatter(eastings, northings, c=np.log10(apparent_stress), edgecolor='k', cmap="nipy_spectral",zorder=22)
        ax2.scatter(c_easting, c_northing, c=np.log10(stress), edgecolor='k', cmap="magma",zorder=22)
        ax3.scatter(c_easting, c_northing, c=np.log10(strain), edgecolor='k', cmap="viridis",zorder=22)   
        ax4.scatter(c_easting, c_northing, c=np.log10(stiffness), edgecolor='k', cmap="nipy_spectral",zorder=22,vmin=4,vmax=16)
        x1,x2 = ax1.get_xlim()
        y1,y2 = ax1.get_ylim()
        for ax in ax1,ax2,ax3,ax4:
            for well in wells.values():
                ax.plot(well['easting'],well['northing'],color='firebrick',lw=3,zorder=11)
                ax.plot(well['easting'],well['northing'],color='0.2',lw=2,zorder=12)
                for well_stage,perf in {k:v for k,v in well.items() if k.isdigit()}.items():
                    stage_eastings = [v for k,v in perf.items() if 'easting' in k]
                    stage_northings = [v for k,v in perf.items() if 'northing' in k]
                    if well_stage==stage.split('_S')[1]:
                        ax.plot(stage_eastings,stage_northings,lw=9,color='k',zorder=14)
                        ax.plot(stage_eastings,stage_northings,lw=7,color='0.8',zorder=15)
                    else:
                        ax.plot(stage_eastings,stage_northings,lw=7,color='k',zorder=15)                      
            ax.set_xlim(x1,x2)
            ax.set_ylim(y1,y2)
        ax1.set_title('events')
        ax2.set_title('clusters\nstress')
        ax3.set_title('clusters\nstrain')
        ax4.set_title('clusters\nstiffness')
        fig.savefig('figures//'+stage+'.png')
        plt.close(fig)