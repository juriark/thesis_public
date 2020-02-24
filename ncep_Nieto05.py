#TODO: download entire_NH data for all variables
#TODO: update step2 and step3

#from netCDF4 import Dataset
import netCDF4 as nc4
import numpy as np

import f_ncep_Nieto05 as myFunc  
import time
import sys
import copy
import pickle
#import cython

R_gas           = 287.058 # [J/(kg*K)] source: wikipedia.org/wiki/Gas_constant (8.11.19, under R_specific for dry air)
G_acc           = 9.80665 # [m/s**2] 
p_heights       = np.array([300, 200])  # has to be decreasing numbers 
years           = np.arange(1958,1998)   # years of data files     #TODO: not np.arange(1958,1999)???
res             = 1     # resolution of grid (in this case: indices are used, not degree. That's why resolution is 1)
event_area      = {}
event_pt        = {}
finished_events = []


# Create netCDF Dataset for writing
step1grp = nc4.Dataset('step1_results.nc','w',format='NETCDF4') #step1grp =; und Zeile danach löschen
step1grp.createDimension('lon',144)
step1grp.createDimension('lat',37)
step1grp.createDimension('time',None)

# Create variables
longitudenc = step1grp.createVariable('Longitude','f4','lon')
latitudenc  = step1grp.createVariable('Latitude','f4','lat')
colnc       = step1grp.createVariable('COL','i4',('time','lat','lon'))
sumcol      = step1grp.createVariable('sumCOL','i4',('lat','lon'))
timenc      = step1grp.createVariable('Time','i4','time')

# Passing data into variables
longitudenc[:]  = np.arange(0,360, 2.5)
latitudenc[:]   = np.arange(90, -2.5, -2.5)

# Adding attributes
step1grp.description   = 'COL points after step 1.2.'

#Add local attributes to variable instances
longitudenc.units   = 'degrees east'
latitudenc.units    = 'degrees north'
timenc.units        = 'hours since Jan 01, 1800'
colnc.units         = 'boolean'
nc_cntr = 0




for year_idx in years:  # files come seperated by years -> open files of current year
    print(year_idx)
    
    geopotential_200    = nc4.Dataset('data_NH/ncepncar_' + str(year_idx) + '_200hPa_hgt.nc','r')
    geopotential_300    = nc4.Dataset('data_NH/ncepncar_' + str(year_idx) + '_300hPa_hgt.nc','r')
    zonal_wind_200      = nc4.Dataset('data_NH/ncepncar_' + str(year_idx) + '_200hPa_uwnd.nc','r')
    temp_200            = nc4.Dataset('data_NH/ncepncar_' + str(year_idx) + '_200hPa_air.nc','r')
    temp_300            = nc4.Dataset('data_NH/ncepncar_' + str(year_idx) + '_300hPa_air.nc','r')
    time                = geopotential_200.variables['time'].size
    
    for timestep in range(time):   #starts with 0, ends with time-1
        #print(' ')
        #print('t = ', timestep, ' __________________________________________')
        
        # get data of current time step
        gp200_timestep      = np.array(geopotential_200.variables['hgt'][timestep,:,:])
        gp300_timestep      = np.array(geopotential_300.variables['hgt'][timestep,:,:])
        zonal_wind_timestep = np.array(zonal_wind_200.variables['uwnd'][timestep,:,:])
        temp_time200        = np.array(temp_200.variables['air'][timestep,:,:])
        temp_time300       = np.array(temp_300.variables['air'][timestep,:,:])
        #temp_profile        = np.array([temp_time200[0,:,:], temp_time300[0,:,:]])
        timestamp           = int(temp_200.variables['time'][timestep])      # hrs since 01-01-1990, timestamp of current timestep
        prev_timestamp      = int(temp_200.variables['time'][timestep-1])    # timestamp for previous timestep
        #print(timestamp)
        
        
        # [step 1.1] 200hPa minimum
        step1_1_boolMask = myFunc.isLocalMin(gp200_timestep, 10, 6)
        print('Number of potential COL points after step[1.1]: ', np.sum(step1_1_boolMask))
        if not step1_1_boolMask.any():
            print('no potential COL points for this time step')
            continue
        
        #TODO:m**2 s**-2?
        
        # [step 1.2] condition of isolated cyclonic vortex
        step1_2_boolMask = myFunc.isIsolatedVortex(step1_1_boolMask, zonal_wind_timestep)
        print('Number of potential COL points after step[1.2]: ', np.sum(step1_2_boolMask))
        if not step1_2_boolMask.any():
            print('no potential COL points for this time step')
            continue
        
        # [step 2] equivalent thickness
        step2_boolMask = myFunc.isEquivalentThickness(step1_2_boolMask, temp_time200[0,:,:], temp_time300[0,:,:], R_gas, G_acc, [200, 1000]) #temp1000
        print('Number of potential COL points after step[2]: ', np.sum(step2_boolMask))
        if not step2_boolMask.any():
            print('no potential COL points for this time step')
            continue
        # TODO woher kommt T_mean? is it really the mean of the profile along a few vertical values? 
        
        # [step 3] Thermal Front Parameter - change of temperature gradient in direction of temperature gradient
        step3_boolMask = myFunc.isTFP(step2_boolMask, temp_time200)
        print('Number of potential COL points after step[3]: ', np.sum(step3_boolMask))
        if not step3_boolMask.any():
            print('no potential COL points for this time step')
            continue
        
        # TODO irgendeinen break einführen, dass er nicht auf nachbarschaft prüft, wenn nur ein potentieller col punkt rauskommt 
        
        # Fill .nc Dataset
        timenc[nc_cntr]     = timestamp
        colnc[nc_cntr,:,:]  = step1_2_boolMask
        nc_cntr             = nc_cntr + 1
        
        
        '''
        # [step 4] Further Restrictions
        col_coordinate = np.array(np.argwhere(step1_2_boolMask == 1)) # get coordinates of COL points #step3_boolMask
        
        
        #print(col_coordinate.size/2,' COL points found at indices: ')
        #print(col_coordinate)
        #print('______________________________________________________')
        #print('Check for neighbourhood:')
        
        
        current_neighs = []
        error = 0
        for i in np.arange(col_coordinate.shape[0]):    #TODO: geht auch 'for i in col_coordinate' oder 'len(col_coordinate)'
            # Get neighbours
            current_pt = col_coordinate[i,:]
            neighbour_idx = myFunc.getNeighbour(np.array([current_pt]), res) # get indices for neighbouring gridpoints
            
            
            #print(' ')
            #print(' ')
            #print('Neighbourhood of ', current_pt, ' is ')
            #print(neighbour_idx)
            #print(' ')        
            
            
            # Check for neighbourhood
            for j in col_coordinate:
                #print('_____________________________________________________')
                #print('col_coordinate to check = ', j)
                
                if (j[0] == current_pt[0]) & (j[1] == current_pt[1]):   #was ist ndim von j und current_pt?
                    #print('col_coordinate is the current_pt -> continue')
                    continue
                elif myFunc.isNeighbour(j, neighbour_idx):
                    current_neighs.append([[j],'is a neighbour of', [current_pt]])
                    #print('col_coordinate ', j, ' is a neighbour')
                elif not myFunc.isNeighbour(j,neighbour_idx):
                    continue
                    #print('col_coordinate ', j, 'is not a neighbour')
                else:
                    #print('ERROR')
                    error = error + 1
        
        # Print final result of current_neighs:
        if bool(current_neighs):    # only enters if current_neighs isn't empty, meaning there are neighbours
            #for a in range(len(current_neighs)):
            #    print(current_neighs[a])
            #print('error = ', error)
            if not error == 0:
                print('Error occured at t = ', timestep)
                sys.exit()
        
        # Group neighbouring points to define area of event
        events_ts = {}
        if not bool(current_neighs):    # Case 1: no neighbours within full_arr -> write full_arr directly as dict
            name = []
            vals = []
            for a in np.arange(len(col_coordinate)):
                name.append('event' + str(a))
                vals.append(np.array([col_coordinate[a,:]]))
            zipObj      = zip(name, vals)       # Create a zip object from two lists (tuples of name and value)
            events_ts   = dict(zipObj)          # Create a dictionary from zip object
        else:                           # Case 2: neighbours -> group neighbours as same events, and then write the dict
            # Convert current_neighs from list to np.array
            current_neighs = np.array(current_neighs)
            neigh_arr = np.zeros([len(current_neighs),4])
            for b in np.arange(len(current_neighs)):    #0 to 8
                neigh_arr[b,0:2] = np.array(current_neighs[b,0])
                neigh_arr[b,2:4] = np.array(current_neighs[b,2])
            neigh_arr = myFunc.removeDuplicates(neigh_arr)
            neigh_arr = neigh_arr.astype(int)
            
            # Do the grouping
            events_ts = myFunc.getEvents(col_coordinate,neigh_arr, False)
        
        # Write event of this timestep either as new events, or as already ongoing event
        if bool(event_area):  # enters if event_area has entries
            sorted_events_output = myFunc.sortIntoEventArea(events_ts, event_area, timestamp, prev_timestamp, finished_events)
            event_area      = copy.deepcopy(sorted_events_output[0])  
            finished_events = copy.deepcopy(sorted_events_output[1])
        else:   # enters if event_area is empty (only at first timestep)
            for layer_event in events_ts.keys():
                event_area[layer_event] = {str(timestamp): events_ts[layer_event]}
        
        output_delOneDayEv = myFunc.deleteOneDayEvents(event_area, finished_events) #code could be more efficient, if it only takes finished_events of current timestep (currently takes finished_events of ALL timesteps). Is even more efficient, if finished_events in sortIntoEventArea is changed accordingly
        event_area = output_delOneDayEv[0]
        finished_events = output_delOneDayEv[1]
        
        # Prepare for next timestep
        del events_ts

    #TODO: check if event_area always deepcopies correctly (has to deepcopy. copy itself is merely some kind of pointer)    
    #TODO: event_point erst ganz zum Schluss erstellen, wenns einmal durch die gesamte Zeitreihe gelaufen ist
    #TODO: southern hemisphere uwnd has to be in other direcntion!
    #TODO: Exclude systems that last only one day


'''
# Close .nc file
col_pts = np.array(step1grp.variables['COL'][:])
sumcol[:,:] = np.array(np.sum(col_pts,axis=0))
step1grp.close()

'''
# open a file, where you want to store the data
file = open('step1_results_entireNH_4ptadjacent', 'wb')
# dump information to that file
pickle.dump(event_area, file)
# close the file
file.close()
'''