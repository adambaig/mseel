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
            wells[well]['northing'][i_line] = float(lspl[12])
            wells[well]['easting'][i_line] = float(lspl[13])
        for stage in [k for k in report.keys() if k.isdigit()]:
            wells[well][stage] = {}
            for cluster,td in {k:float(v) for k,v in report[stage].items() if 'TD' in k and v!=''}.items():
                easting = np.interp(td,wells[well]['md'], wells[well]['easting'])
                northing = np.interp(td,wells[well]['md'], wells[well]['northing'])
                elevation= np.interp(td,wells[well]['md'], wells[well]['northing'])                
                wells[well][stage][cluster+' easting'] = easting
                wells[well][stage][cluster+' northing'] = northing
                wells[well][stage][cluster+' elevation'] = elevation
            
            
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
    
 

    
    
    