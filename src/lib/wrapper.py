#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 13:34:59 2019

@author: davidmorante
"""
#
# READ OEM FILE
#
import pandas as pd
import numpy as np
import json


def oem_wrapper(input_file):
    import pandas as pd
    import numpy as np
    import json

    with open(input_file) as f:
            pos = 0
            cur_line = f.readline()
            while not cur_line.startswith("META_STOP"):
                pos = f.tell()
                cur_line = f.readline()
            cur_line = f.readline()
            pos = f.tell()
            f.seek(pos)
            data2 = np.array(pd.read_csv(f, header=None, sep=' '))
    
    f.close()      
    time = data2[0:,0] + data2[0:,1]/(3600*24)
    R    = data2[0:,2:5] 
    V    = data2[0:,5:9] 
    
    return time, R, V 

def modes_wrapper(input_file):
    with open(input_file) as f:
            pos = 0
            cur_line = f.readline()
            while not cur_line.startswith("META_STOP"):
                pos = f.tell()
                cur_line = f.readline()
            cur_line = f.readline()
            pos = f.tell()
            f.seek(pos)
            data2 = np.array(pd.read_csv(f, header=None, sep=' '))
    
    f.close()      
    time = data2[0:,0] + data2[0:,1]/(3600*24)
    #caca  = data2[0,2]
    return time

def json_wrapper(input_file):
    #################################################
    # JSON FILE WRAPPER FOR THE ACTIVITY PROFILE
    #################################################
    with open(input_file) as f:
        data = json.load(f)
    
    f.close()
    level0 = data['children']
    objs = []
    
    for i in range(0,len(level0)):
        objs.append(level0[i]['name'])
        
    #################################################
    # Obtain Spacecraft Orbital Elements
    #################################################
        
    ind   = objs.index('Orbitography') 
    orbitography = level0[ind]['children']
    
    objs2 = []
    for i in range(0,len(orbitography)):
        objs2.append(orbitography[i]['name'])
        
    ind2   = objs2.index('Satellite_orbit') 
       
    orbit = orbitography[ind2]['modesChildren'][0]['values']
    
    objs3 = []
    for i in range(0,len(orbit)):
        objs3.append(orbit[i]['name'])
        
    #
    # RAAN
    #
    ind_raan = objs3.index('RAAN') 
    raan     = float(orbit[ind_raan]['value']) 
    #
    # argument of perigee
    #
    ind_w = objs3.index('argument_of_perigee') 
    w      = float(orbit[ind_w]['value']) 
    #
    # init_cond_on_true_anomaly
    #
    ind_ta = objs3.index('init_cond_on_true_anomaly') 
    ta     = float(orbit[ind_ta]['value']) 
    #
    # inclination
    #
    ind_incl = objs3.index('inclination') 
    incl     = float(orbit[ind_incl]['value']) 
    #
    # Eccentricity (Assumed to 0)
    #
    ind_ecc = objs3.index('eccentricity') 
    ecc= float(orbit[ind_ecc]['value']) 
    #
    # Semimajor axis, if not available search for altitude
    #
    try:
        ind_a = objs3.index('semi_major_axis') 
        a     = float(orbit[ind_a]['value'])
    except:
        ind_alt = objs3.index('altitude') 
        alt     = float(orbit[ind_alt]['value'])
        a     = alt/1000 + 6378; # meters
    #
    # Mission_time
    #
    ind_mission_time = objs3.index('mission_time') 
    mission_time     = float(orbit[ind_mission_time]['value']) 
    mission_time     = mission_time * 2 * 3.1416*(a*a*a/398600.440)**(0.5)
    #
    # Initial Date
    #
    ind_init_date = objs3.index('init_date') 
    init_date     = float(orbit[ind_init_date]['value']) 
    #
    ## Convert to a data Frame
    #
    coe = {
    "a" : [a],
    "e" : [ecc],               
    "incl" : [incl*np.pi/180],  
    "w" : [w*np.pi/180],
    "RA": [raan*np.pi/180],
    "TA": [ta*np.pi/180]
    }
    coe = pd.DataFrame(coe) 
    
    #################################################
    # Obtain Spacecraft Modes
    #################################################
    #  
    ind = objs.index('Satellite') 
    #  
    satellite = level0[ind]['modesChildren']
    #
    mode_name = []
    mode_power = []
    for i in range(0,len(satellite)):
        mode_name.append(satellite[i]['name'])
        mode_power.append(0)
    #
    ## Convert to a data Frame
    #
    scmodes = {
    "Name"  : mode_name,
    "Power" : mode_power
    }
    scmodes = pd.DataFrame(scmodes) 
    #
    #################################################
    # Obtain Ground Stations
    #################################################
    #
    ind = objs.index('Ground_Segment') 
    #  
    ground_segment = level0[ind]['modesChildren']
    #
    ground_station = []
    for i in range(0,len(ground_segment)):
        ground_station.append(ground_segment[i]['name'])
    #
    ground_station_name      = []
    ground_station_latitude  = []
    ground_station_longitude = []
    ground_station_altitude  = []
    ground_station_elevation = []
    #
    for i in range(0,len(ground_station)):
        obj4 = []
        for j in range(0,len(ground_segment[i]['values'])): 
               obj4.append(ground_segment[i]['values'][j]['name'])
               
        ground_station_latitude.append(float(ground_segment[i]['values'][obj4.index('Latitude')]['value'])*np.pi/180)
        ground_station_longitude.append(float(ground_segment[i]['values'][obj4.index('Longitude')]['value'])*np.pi/180)
        ground_station_altitude.append(float(ground_segment[i]['values'][obj4.index('Altitude')]['value']))
        ground_station_elevation.append(float(ground_segment[i]['values'][obj4.index('Elevation')]['value'])*np.pi/180)
        ground_station_name.append(ground_segment[i]['values'][obj4.index('Name')]['value']) 
    #
    ## Convert to a data Frame
    #
    gs = {
    "Latitude" : ground_station_latitude,
    "Longitude" : ground_station_longitude,
    "Altitude" : ground_station_altitude,
    "Elevation" : ground_station_elevation,
    "Name": ground_station_name
    }
    gs = pd.DataFrame(gs)     
            
    return coe, gs , scmodes, init_date, mission_time



