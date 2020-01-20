from datetime import datetime, timedelta
import logging
import numpy as np
from read_inputs import read_and_rotate_stage


def calc_distance_sq(point1,point2):
    return (point1['easting'] - point2['easting'])**2 + (point1['northing'] - point2['northing'])**2 + (point2['elevation'] - point2['elevation'])**2


def nearest_neighbors(events,point,t_win=None,distance_range=None, min_events=None):
    if t_win==None and distance_range==None and min_events==None:
        logging.error('specify t_win, distance_range, and/or min_events')
        return
    if distance_range==None:
        distance_range = 1.e10
    distance_range_sq = distance_range*distance_range
    if len(t_win)==2:
        culled_events = {k:v for k,v in events.items() if 
            (v['time_zero'] > t_win[0] and v['time_zero'] < t_win[1])}
    else:
        culled_events = events
    if range is not None:  # do a first selection before calculating distances
        culled_events = {k:v for k,v in culled_events.items() if 
             (v['easting'] > point['easting'] - distance_range and
              v['easting'] < point['easting'] + distance_range and
              v['northing'] > point['northing'] - distance_range and
              v['northing'] < point['northing'] + distance_range and                         
              v['elevation'] > point['elevation'] - distance_range and
              v['elevation'] < point['elevation'] + distance_range)}
    distance_sq = {}
    for i_event,(event_id,event) in enumerate(culled_events.items()):
        distance_sq[event_id] = calc_distance_sq(event,point)
    distance_sq_array = [(v,k) for k,v in distance_sq.items() if v < distance_range_sq]
    return {event_id[1]:culled_events[event_id[1]] for event_id in distance_sq_array[:min_events]}


def cluster(stage,well,min_events=10,t_win=None,range=50):
    stage_events = read_and_rotate_stage('20','3H')

    fig,ax = plt.subplots()

    t_zero = np.array([timedelta(0,0,v['t0']*1.e6) + 
                       datetime.strptime(v['time'],'%d:%m:%y:%H:%M:%S') 
                       for k,v in stage_events.items() if v['MS_EVENT_TYPE']==0])
    i_sort = np.argsort(t_zero)
    easting = np.array([v['easting'] for k,v in stage_events.items() if v['MS_EVENT_TYPE']==0])[i_sort]
    northing = np.array([v['northing'] for k,v in stage_events.items() if v['MS_EVENT_TYPE']==0])[i_sort]
    h1 = np.array([v['h1'] for k,v in stage_events.items() if v['MS_EVENT_TYPE']==0])[i_sort]
    h2 = np.array([v['h2'] for k,v in stage_events.items() if v['MS_EVENT_TYPE']==0])[i_sort]
    elevation = np.array([v['elevation'] for k,v in stage_events.items() if v['MS_EVENT_TYPE']==0 ])[i_sort]
    
    
    

    

