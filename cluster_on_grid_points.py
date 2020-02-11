import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
from scipy.ndimage.filters import gaussian_filter
from scipy.spatial import ConvexHull, Delaunay
from sklearn.neighbors import  BallTree
from read_inputs import read_rdvs, read_wells

k_neighbours = 5
grid_spacing= 30
events = read_rdvs(perfs=False)

wells= read_wells()

complete_events =  {k:v for k,v in events.items() if v["SP_MAGNITUDE"]>-2}

eastings = [v['easting'] for v in complete_events.values()]
northings = [v['northing'] for v in complete_events.values()]
elevations = [v['elevation'] for v in complete_events.values()]
apparent_stress = [v["SP_ENERGY"]/v["SP_POTENCY"] for v in complete_events.values()]
event_keys= list(complete_events.keys())


east_grid = np.arange(min(eastings), max(eastings), grid_spacing)
north_grid = np.arange(min(northings), max(northings), grid_spacing)
elev_grid = np.arange(min(elevations), max(elevations), grid_spacing)
  

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
neighbour_tree = BallTree(xyz)

ref_elevation= -6380
east_grid = np.arange(min(eastings),max(eastings),grid_spacing)
north_grid = np.arange(min(northings), max(northings), grid_spacing)
elevation_grid = np.arange(min(elevations), max(elevations), grid_spacing)

i_slice = np.where(elevation_grid<ref_elevation)[0][-1]


grid_points = []
for i_east,e_grid in enumerate(east_grid):
    for i_north, n_grid in enumerate(north_grid):
        grid_points.append([e_grid, n_grid, ref_elevation])

ne_grid, nn_grid,nz_grid = len(east_grid), len(north_grid),len(elevation_grid)
mask = np.logical_not(in_hull(grid_points, xyz)).reshape([ne_grid, nn_grid])



stress_map= np.zeros([ne_grid, nn_grid])
strain_map = np.zeros([ne_grid, nn_grid])
stiffness_map = np.zeros([ne_grid, nn_grid])
for i_east,e_grid in enumerate(east_grid):
    for i_north, n_grid in enumerate(north_grid):
        dist, ind = neighbour_tree.query([[e_grid, n_grid, ref_elevation]], k_neighbours)
    #@    if in_hull([e_grid, n_grid,ref_elevation], xyz[ind[0]]):
        volume = ConvexHull(xyz[ind[0]]).volume*3.28**3
        strain_map[i_east, i_north] = sum([complete_events[event_keys[i]]['SP_POTENCY'] for i in ind[0]])/2/volume
        stress_map[i_east, i_north] = sum([complete_events[event_keys[i]]['SP_ENERGY'] for i in ind[0]])/strain_map[i_east, i_north]
        stiffness_map[i_east, i_north] = stress_map[i_east, i_north]/strain_map[i_east, i_north]


fig1, ax1 = plt.subplots(1,figsize=[10,8])
fig2, ax2 = plt.subplots(1,figsize=[10,8])
fig3, ax3 = plt.subplots(1,figsize=[10,8])
fig4, ax4 = plt.subplots(1,figsize=[10,8])
for ax in [ax1,ax2,ax3,ax4]:
    ax.set_aspect('equal')
    ax.set_facecolor('0.95')
strain_masked=ma.masked_array(gaussian_filter(strain_map, 0.8), mask=mask)
stress_masked=ma.masked_array(gaussian_filter(stress_map, 0.8), mask=mask)
stiffness_masked=ma.masked_array(gaussian_filter(stiffness_map,0.8), mask=mask)
apparent_stress_plot = ax1.scatter(eastings, northings, c=np.log10(apparent_stress), edgecolor='k', cmap="nipy_spectral",zorder=2)
strain_plot = ax2.pcolor(east_grid, north_grid, np.log10(strain_masked.T), cmap="nipy_spectral",zorder=2)
stress_plot = ax3.pcolor(east_grid, north_grid, np.log10(stress_masked.T),  cmap="nipy_spectral",zorder=2)   
compliance_plot = ax4.pcolor(east_grid, north_grid, np.log10(stiffness_masked.T), cmap="nipy_spectral",zorder=2)
x1,x2 = ax1.get_xlim()
y1,y2 = ax1.get_ylim()
for ax in ax1,ax2,ax3,ax4:
    for well in wells.values():
#        ax.plot(well['easting'],well['northing'],color='firebrick',lw=3,zorder=11)
        ax.plot(well['easting'],well['northing'],color='0.2',lw=2,zorder=12,alpha=0.5)
        for well_stage,perf in {k:v for k,v in well.items() if k.isdigit()}.items():
            stage_eastings = [v for k,v in perf.items() if 'easting' in k]
            stage_northings = [v for k,v in perf.items() if 'northing' in k]
            ax.plot(stage_eastings,stage_northings,lw=7,color='k',zorder=15,alpha=0.5)                      
    ax.set_xlim(x1,x2)
    ax.set_ylim(y1,y2)
for ax in [ax1,ax2,ax3,ax4]:
    ax.set_xlabel('easting (ft)')
    ax.set_ylabel('northing (ft)')
cb1 = fig1.colorbar(apparent_stress_plot)
cb2 = fig2.colorbar(strain_plot)
cb3 = fig3.colorbar(stress_plot)
cb4 = fig4.colorbar(compliance_plot)
ax1.set_title('events')
ax2.set_title('clusters\nstrain')
ax3.set_title('clusters\nstress')
ax4.set_title('clusters\nstiffness')
cb1.set_label('log apparent stress (Pa)')
cb2.set_label('log strain')
cb3.set_label('log stress (Pa)')
cb4.set_label('log stiffness (Pa)')