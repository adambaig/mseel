import glob
import numpy as np


def read_rdvs(events=True, perfs=True):
    dr = 'C:\\Users\\Adam\\Documents\\mseel\\analysis\\?H\\RDV\\'
    if (events and perfs):
        ms_files = glob.glob(dr+'*')
    elif (events):
        ms_files = glob.glob(dr+'*S??F*')    
    elif (perfs):
        ms_files = glob.glob(dr+'*S??P*')
    else:
        return          
    events = {}
    coordinates = {'X': 'northing', 'Y': 'easting', 'Z': 'elevation', '0': 't0'}
    for ms_file in ms_files:
        f = open(ms_file)
        head= f.readline()
        lines = f.readlines()
        f.close()
        well = ms_file.split('_')[-3]
        stage = ms_file.split('_')[-2].split('-FINAL')[0]
        parameters = head.split('#NULL')[0].split()[1:]
        well_stage = well+'_'+stage
        for line in lines[2:]:
            lspl = line.split()
            event_id = well_stage+'_'+lspl[1].zfill(6)
            events[event_id] = {}
            for i_parameter, parameter in enumerate(parameters):
                if i_parameter==0:
                    events[event_id]['time'] = lspl[0]
                elif parameter in  ['MS_EVENT_TYPE', 'SP_COMPONENTS_USED', 'SP_RECEIVERS_USED']:
                    events[event_id][parameter] = int(lspl[i_parameter])
                elif 'QC_LOC_' in parameter:
                    coord = coordinates[parameter[-1]]
                    events[event_id][coord] = float(lspl[i_parameter])
                elif i_parameter!=1:
                    events[event_id][parameter] = float(lspl[i_parameter])                               
    return(events)
       
def read_wells():
    wells = {}
    for well_csv in glob.glob('C://Users//Adam//Documents//mseel//analysis//Wells//*.csv'):
        well= well_csv.split('\\')[-1].split('.')[0]
        wells[well] = {}
        f = open(well_csv)
        for ii in range(22):
            dum = f.readline()
        lines = f.readlines()
        f.close()
        n_lines = len(lines)
        report = read_stage_reports(well)
        wells[well] = {
                'md': np.zeros(n_lines),
                'elevation': np.zeros(n_lines),
                'northing': np.zeros(n_lines),
                'easting': np.zeros(n_lines),                
                }
        for i_line,line in enumerate(lines):
            lspl = line.split(',')
            wells[well]['md'][i_line] = float(lspl[1])
            wells[well]['elevation'][i_line] = -float(lspl[5])
            wells[well]['northing'][i_line] = float(lspl[13])
            wells[well]['easting'][i_line] = float(lspl[12])
        for stage in [k for k in report.keys() if k.isdigit()]:
            wells[well][stage] = {}
            for cluster,td in {k:float(v) for k,v in report[stage].items() if 'TD' in k and v!=''}.items():
                easting = np.interp(td,wells[well]['md'], wells[well]['easting'])
                northing = np.interp(td,wells[well]['md'], wells[well]['northing'])
                elevation= np.interp(td,wells[well]['md'], wells[well]['northing'])
                position = cluster.replace('Top','Heel').replace('Bottom','Toe')
                wells[well][stage][position+' easting'] = easting
                wells[well][stage][position+' northing'] = northing
                wells[well][stage][position+' elevation'] = elevation
    return(wells)
        
def read_stage_reports(well):
    f = open('C:\\Users\\Adam\\Documents\\mseel\\analysis\\'+well+'\\NNE MIP '+
             well+' Final Frac Summary.txt')
    for ii in range(14):
        dum = f.readline()
    line_stage = f.readline()
    dum = f.readline()
    report = {}
    stage_str = []
    for stage in line_stage.split('\t')[2:]:
        report[stage.replace('\n','')] = {}
        stage_str.append(stage.replace('\n',''))
    for line in f.readlines():
        lspl = line.split('\t')
        attribute = lspl[0]
        for i_value, value in enumerate(lspl[2:]):
            report[stage_str[i_value]][attribute] = value.replace(',','').replace('\n','').replace('\"','')
    f.close()
    return(report)
    
wells = read_wells()
events = read_rdvs()
    
def read_and_rotate_stage(stage,well,azimuth = None):
    perfs = wells[well][stage]
    stage_events = {k:v for k,v in events.items() if (well+'_S'+stage) in k}
    perf_eastings = np.array([v for k,v in perfs.items() if 'easting' in k])
    perf_northings = np.array([v for k,v in perfs.items() if 'northing' in k])
    perf_elevations = np.array([v for k,v in perfs.items() if 'elevation' in k])
    avg_easting = np.average(perf_eastings)
    avg_northing = np.average(perf_northings)
    avg_elevation = np.average(perf_elevations)
    if azimuth == None:
        azimuth = np.arctan2(perf_northings[-1] - perf_northings[0],
                             perf_eastings[-1] - perf_eastings[0])
    else:
        azimuth=0
    R = np.array([[np.cos(azimuth), -np.sin(azimuth)],
                  [np.sin(azimuth), np.cos(azimuth)]])
    for event in stage_events:
        rel_east = stage_events[event]['easting'] - avg_easting
        rel_north = stage_events[event]['northing'] - avg_northing
        rel_elevation = stage_events[event]['elevation'] - avg_elevation
        h1,h2 = R@np.array([rel_east,rel_north]).T
        stage_events[event]['h1'] = h1
        stage_events[event]['h2'] = h2
        stage_events[event]['rel_elevation'] = rel_elevation
    return stage_events

def read_treatment_rdvs(well): 
    return
    