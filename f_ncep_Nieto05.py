import numpy as np
import math
import sys
import copy

#####################################################################################################
# Print only. Gives first glimpse of nc data (variables, dimensions)
#####################################################################################################
def exploreData(nc_file):
    print('start#########################################')
    print('netCDF file format: ' + nc_file.file_format) # displays netCDF file format
    print('\n')

    print('netCDF variables: ', nc_file.variables.keys()) # view netCDF file's variables
    for i in nc_file.variables: # displays variables, current shape and units
        print(i, ': ', nc_file.variables[i].shape, ', ', nc_file.variables[i].units)
    print('\n')

    print('netCDF dimensions:')
    print(nc_file.dimensions.keys()) # python dictionary calls allow viewing of the netCDF's file dimensions
    print('##########################################end')


#####################################################################################################
# Finds the local minima of 3x3 subgrid on a global surface
#
# Input:    myArray             - numpy array of geopotential field
#           threshold           - geopotential difference with surrounding points, that needs to be 
#                                   exceeded
#           numel_neighbours    - number of neighbours that have to be smaller than central point in 
#                                   order for it to be a local minimum
# Output:   boolOut             - bool-mask of either 0 (not a potential COL center) or 1 (potential 
#                                   COL center)
#####################################################################################################
def isLocalMin(myArray, threshold, numel_neighbours):
    myArray = myArray[0,:,:]
    dim1 = myArray.shape
    lat_ini = np.transpose(np.array([myArray[:,0]]))
    lat_end = np.transpose(np.array([myArray[:,dim1[1]-1]]))
    expandArray = np.append(lat_end, np.append(myArray, lat_ini, axis=1), axis=1) # append values at lat=0 and lat=357.5 to ensure evaluation over lat=0
    dim2 = expandArray.shape
    boolOut = np.zeros(dim1) # initialize output array
    neighbour = np.zeros(8)       
    #for lat in np.arange(1,dim2[0]-1):                   #beginnt wirlkich ab 1? wird dann nicht lat[0] gar nicht ausgewertet?
    for lat in np.arange(8,29):  #idx=8 entspricht bei 2.5°res Breitengrad 70 (sind 20° Abstand von 90°N)
        for lon in range(1, dim2[1]-1):    
            centerPoint = expandArray[lat,lon]
                
            # get neighbours
            neighbour[0] = expandArray[lat-1  ,lon-1]
            neighbour[1] = expandArray[lat-1  ,lon  ]
            neighbour[2] = expandArray[lat-1  ,lon+1]
            neighbour[3] = expandArray[lat    ,lon-1]
            neighbour[4] = expandArray[lat    ,lon+1]
            neighbour[5] = expandArray[lat+1  ,lon-1]
            neighbour[6] = expandArray[lat+1  ,lon  ]
            neighbour[7] = expandArray[lat+1  ,lon+1]
            
            #if lon == dim[1]:     # necessary to connect 359.75° and 0° longitude
            #    neighbour[2] = expandArray[lat-1  ,0]
            #    neighbour[4] = expandArray[lat    ,0]
            #    neighbour[7] = expandArray[lat+1  ,0]
            #else:
            #    neighbour[2] = expandArray[lat-1  ,lon+1]
            #    neighbour[4] = expandArray[lat    ,lon+1]
            #    neighbour[7] = expandArray[lat+1  ,lon+1]

            #check for minimum
            diff = neighbour[:] - centerPoint
            cntr = 0
            for n in range(8):
                if diff[n] >= threshold:
                    cntr += 1
            if cntr >= numel_neighbours:
                boolOut[lat,lon-1] = 1    #here lon-1???              
    return(boolOut)


#####################################################################################################
# Checks set of potential COL center for condition of isolated cyclonic vortex. Here, this is means 
# a change of zonal wind direction at at least one of the two northern grid points.
#
# Input:    myBoolMask          - bool-mask of potential COL centers (as in Output of 
#                                   def findMinima())
#           zonalWindArray      - numpy array of zonal wind speed u (only sign is of interest)
# Output:   myBoolMask          - bool-mask of potential COL points fulfilling isolated cyclonic 
#                                   vortex condition
#####################################################################################################
def isIsolatedVortex(myBoolMask, zonalWindArray):
    zonalWindArray = zonalWindArray[0,:,:]
    ones_idx = np.argwhere(myBoolMask == 1) # get indices of potential COL points
    for i in ones_idx:
        lat = i[0]
        lon = i[1]
        my_u        = np.sign(zonalWindArray[lat, lon]) # get signs
        north_u1    = np.sign(zonalWindArray[lat-1, lon])
        north_u2    = np.sign(zonalWindArray[lat-2, lon])
        if (my_u == north_u1 * -1) or (my_u == north_u2 * -1):  # check for condition
            continue
        else: # if isolated vortex condition not fulfilled, replace 1 with 0
            myBoolMask[lat,lon] = 0
    return(myBoolMask)


####################################################################################################
# Checks for condition of heigher equivalnet thickness eastward of the potential COL central point
# Equation: h = R*mean(T_v)/g * ln(p_1/p_2)  
#
# Input:    myBoolMask          - bool-mask of potential COL centers (as in Output of 
#                                   def isIsolatedVortex())
#           temp200             - temperature at 200 hPa
#           temp1000            - temperature at 1000 hPa
#           R_gas               - specific gas constant for dry gas
#           g_acc               - gravitational acceleration
#           p_heights           - list containing pressure surface heights of heighest and lowest 
#                                   level
# Output:   myBoolMask          - bool-mask of potential COL points fulfilling equivalent thickness 
#                                   condition
#####################################################################################################
def isEquivalentThickness(myBoolMask, temp200, temp1000, R_gas, g_acc, p_heights):
    ones_idx = np.argwhere(myBoolMask == 1) # get indices of potential COL central points
    for i in ones_idx:
        lat     = i[0]
        lon     = i[1]
        eastlon = i[1]+1 if not (lon == temp200.shape[1]-1) else 0  # else-statement, to account for discontinuity 
        
        t_virt_central      = temp200[lat,lon] * (1 + 0.61)
        t_mean_central      = 
        t_mean_east         =       
        thickness_central   = R_gas * t_mean_central / g_acc * np.log(p_heights[0]/p_heights[1])
        thickness_east      = R_gas * t_mean_east    / g_acc * np.log(p_heights[0]/p_heights[1])
        if thickness_central >= thickness_east:
            myBoolMask[lat,lon] = 0
        else:
            continue
    return(myBoolMask)        


#####################################################################################################
# Checks for condition of Thermal Front Parameter (TFP). TFP eastward of central point must be      
# higher.                                                                                           
# TFP = -grad |grad(t)| * ( grad(T) / |grad(T)| )                                                       
#       ----left term---   ------right term------                                                   
#                                                                                                   
# Input:    myBoolMask          - bool-mask of potential COL centers                                
#                                   (as in Output def isEquivalentThickness())                      
#           tempArray           - temperature numpy array                                           
# Output:   myBoolMask          - bool-mask of potential COL points fulfilling TFP condition        
#####################################################################################################
def isTFP(myBoolMask, tempArray):
    tempArray = tempArray[0,:,:]
    grad_T = np.array(np.gradient(tempArray))
    abs_grad_T = np.zeros(tempArray.shape)  # initialize array
    for lat in range(tempArray.shape[0]):
        for lon in range(tempArray.shape[1]):
            abs_grad_T[lat,lon] = math.sqrt(grad_T[0][lat,lon]**2 + grad_T[1][lat,lon]**2)   # length of vector grad(T) at every gridpoint

    left_term = np.gradient(abs_grad_T)
    ones_idx = np.argwhere(myBoolMask == 1) # get indices of potential COL central points  

    for i in ones_idx:
        lat     = i[0]
        lon     = i[1]
        eastlon = i[1]+1 if not (lon == myBoolMask.shape[1]-1) else 0
        
        grad_T_cent = np.array([grad_T[0][lat, lon], grad_T[1][lat,lon]])   # [2x1 vector]
        grad_T_east = np.array([grad_T[0][lat, eastlon], grad_T[1][lat, eastlon]])

        abs_grad_T_cent = math.sqrt(grad_T_cent[0]**2 + grad_T_cent[1]**2) # scalar
        abs_grad_T_east = math.sqrt(grad_T_east[0]**2 + grad_T_east[1]**2)

        right_term_cent = grad_T_cent / abs_grad_T_cent    #[2x1] vector
        right_term_east = grad_T_east / abs_grad_T_east

        left_term_cent = np.array([[left_term[0][lat, lon], left_term[1][lat, lon]]])
        left_term_east = np.array([left_term[0][lat, eastlon], left_term[1][lat, eastlon]])

        tfp_cent = np.dot(left_term_cent,right_term_cent) # scalar
        tfp_east = np.dot(left_term_east,right_term_east)

        if tfp_cent > tfp_east:
            myBoolMask[lat,lon] = 0
        else:
            continue
    return(myBoolMask)


#####################################################################################################
# Gets the neighbourhood coordinates of a matrix
#
# Input:    myArray             - numpy 1d array to create neighbourhood for. If input-array 
#                                   contains multiple grid points, these must be neighbours already.
#                                   Input shape is (2*number of gridpoints,)
#           res                 - resolution of grid (even spacing between grid points assumed)
# Output:   myNeighb            - [mx2] numpy array containing coordinates of neighbours (excluding
#                                   myArray-gridpoints and duplicates)
#####################################################################################################
def getNeighbour(myArray, res):
    numel_points = myArray.shape[0]    
    neighbours = np.zeros([8*numel_points,2])

    for i in range(numel_points):   #compute neighbours for gridpoints seperately
        col_point = myArray[i]
        neigh_timestep = np.array([   
            [col_point[0]-res, col_point[1]-res ],
            [col_point[0]-res, col_point[1]     ],
            [col_point[0]-res, col_point[1]+res ],
            [col_point[0]    , col_point[1]-res ],
            [col_point[0]    , col_point[1]+res ],
            [col_point[0]+res, col_point[1]-res ],
            [col_point[0]+res, col_point[1]     ],
            [col_point[0]+res, col_point[1]+res ]])
        neighbours[8*i : 8*(i+1), :] = neigh_timestep

    myNeigh = np.array(np.unique(neighbours, axis=0)) # delete duplicates
    for i in range(numel_points):   # delete myArray-gridpoints
        xarr = myArray[i,0]
        yarr = myArray[i,1]
        xngh = myNeigh[:,0]
        yngh = myNeigh[:,1]
        idx  = np.argwhere( (xngh == xarr) & (yngh == yarr) )
        myNeigh = np.delete(myNeigh, idx, 0)
    return myNeigh    


#####################################################################################################
# Checks for equal values within myArray and myNeighbours-array. If there is, function returns True.
#
# Input:    myArray             - numpy Array containing x and y coordinate(s) 
#           myNeighbours        - numpy Array containing neighbourhood
# Output:   boolOut             - returns True if neighbours, and False if not
#####################################################################################################
def isNeighbour(myArray, myNeighbours):
    #TODO wie handelt man die Nachbarschaft über Breitengrad 0 hinweg?
    if myArray.ndim == 1:
        myArray = np.array([myArray])   # transforms 1d in 2d array
    
    # get myArray in shape (n,2) - Vektorschreibweise
    if (not myArray.shape[1] == 2) & (myArray.shape[0] == 2):
        myArray = np.transpose(myArray)
    elif (not myArray.shape[1] == 2) & (not myArray.shape[0] == 2):
        print('Error with myArray.shape. Expected shape (n,2) or (2,n), but found shape ', myArray.shape)
        sys.exit()
    
    # checks if myArray[0] is contained within myNeighbours
    for i in np.arange(len(myNeighbours)):
        if (myArray[0,0] == myNeighbours[i,0]) & (myArray[0,1] == myNeighbours[i,1]):
            return True #TODO if it reaches return, does it really exit the function?
        else:
            continue
    
    if myArray.shape[0] == 1:
        return False
    else:
        isNeighbour(myArray[1:myArray.shape[0]], myNeighbours) # recursive function call


#####################################################################################################
# Remove duplicate values in neighbourhood matrix
#
# Input:    neigh_arr           - [n,4] numpy array. Columns give neighbouring grid points
# Output:   singles             - same as input, but duplicate coordinate-pairs are removed
#####################################################################################################
def removeDuplicates(neigh_arr):
    duplicates  = int(len(neigh_arr)/2)
    #print(duplicates)
    mirrored    = np.zeros(neigh_arr.shape)
    mirrored[0:duplicates] = neigh_arr[0:duplicates]
    for i in np.arange(duplicates, len(neigh_arr)):
        mirrored[i,0:2] = neigh_arr[i,2:4]
        mirrored[i,2:4] = neigh_arr[i,0:2]
    singles = np.unique(mirrored, axis=0)
    return(singles)


#####################################################################################################
# Writes pairs of neighbours into an array, where every row contains entire area of one COL event.
# A treshhold of -10000 is used to minimize redundant function call.
#
# Input:    full_arr            - [n,2] numpy array. Contains the full array of COL points
#           neigh_arr           - [n,2] numpy array. Contains information on neighbourhood. Lines
#                                 give adjacent points
#           recursion           - Always set to False when calling the function from outside the
#                                 function. Only True within the function itself.
# Output:   events4timestep     - dict. Keys give the different events for the current timestep. 
#                                 Values are [n,2] numpy arrays, with x and y in the columns.
#####################################################################################################
def getEvents(full_arr, neigh_arr, recursion):
    event_cntr      = 0
    tresh           = -10000
    events4timestep = {}
    neigh_len       = len(neigh_arr)
    
    if neigh_arr.shape[1] != 2:                                                                         # making sure neigh_arr is in the necessary shape
        neigh_vec = np.append(neigh_arr[:,0:2], neigh_arr[:,2:4], axis=0)
    else:
        neigh_vec = neigh_arr
    
    # loop through the full array of COL points
    for x in np.arange(len(full_arr)):
        if np.array_equal(full_arr[x], np.array([tresh,tresh])):                                        # continue, if current val = threshold
            continue
        current_event = np.array([full_arr[x]])
        
        for i in np.arange(len(neigh_vec)):                                                             # loop to check if full_arr has a neighbour
            if np.array_equal(neigh_vec[i], current_event[0]):                                          # full_arr_val is in neigh_arr, which means it has neighbour
                neigh_idx       = i - neigh_len                                                         # gets index of neighbour
                current_event   = np.append(current_event, [neigh_vec[neigh_idx]],axis=0)               # appends the neighbour to the current_event var
                
                # set this neighbour in neigh_arr and neigh_vec to tresh (don't delete them, or loop-indices need to be altered as well)
                neigh_arr[neigh_idx]    = tresh                                                         # update neigh_arr by eliminating 'used' neighbourhood information (to get rid of redundancy)
                neigh_vec               = np.append(neigh_arr[:,0:2], neigh_arr[:,2:4], axis=0)
                
                # recursive function call with the newly identified neighbour
                recursion_list  = getEvents([current_event[len(current_event)-1]], neigh_arr, True) 
                current_event   = np.append(current_event, recursion_list[0], axis=0)
                neigh_arr       = recursion_list[1]
                neigh_vec       = np.append(neigh_arr[:,0:2], neigh_arr[:,2:4], axis=0)
                
        if recursion:
            continue
        else:   
            name = 'event' + str(event_cntr)
            event_cntr = event_cntr + 1
            events4timestep[name] = np.unique(current_event, axis=0)
            
            # delete the as neighbour identified val from full_arr (is already handled within recursive function call)
            for i in np.arange(len(full_arr)):
                for j in np.arange(len(current_event)):
                    if np.array_equal(full_arr[i], current_event[j]):
                        full_arr[i] = tresh
                        
    # return statement
    if recursion:
        return [current_event, neigh_arr]
    else:
        return events4timestep #TODO: np array in events4timestep are shape (n,). make 2d arrays!


#####################################################################################################
# Takes an event dictionary of one timestep, and sorts the data into an event dict, containing
# events over multiple areas and timesteps.
# event_area is in the shape    event_area = {
#                                   'event1': {
#                                           't_n': np.array of coordinates of a COL event at t_n,
#                                           't_n+1': np.array of coordinates of a COLevent at t_n+1,
#                                               ...
#                                           }
#                                   'event2': {
#                                           't_m': np.array of coordinates of a COL event at t_m,
#                                           't_m+1': np.array of coordinates of a COLevent at t_m+1,
#                                               ...
#                                           }
#                                       ...
#                                       }
#
# Input:    event_ts                - dict, events of current timestep
#           event_area              - dict, events of all previous timesteps
#           current_ts              - int, current timestamp, hours since 01-01-1990
#           prev_ts                 - int, timestamp of previous timestep, hours since 01-01-1990
#           finished_events         - list, containing the names of past events, that have ended already
# Output:   event_area              - updated event_area dict
#           finished_events_update  - updated finished_events list
#####################################################################################################
def sortIntoEventArea(event_ts, event_area, current_ts, prev_ts, finished_events):
    updated_event_area = copy.deepcopy(event_area)
    cntr = 1
    keys = [keys for keys in event_area.keys() if keys not in finished_events]
    
    for name_curr in event_ts.keys():   # loop trough events of current timestep
        current_event       = event_ts.get(name_curr)     # np array
        notBeenAttributed   = True
        
        #for name_prev in event_area.keys():    # loop trough events found in previous timesteps  
        for name_prev in keys:
            existingEvent   = False
            
            #if name_prev in finished_events:    # Check if name_prev is an already finished event. If yes, continue
            #    continue        #TODO: code could be made more efficent if skips finished events already before
            
            prev_event      = event_area.get(name_prev)  # pair of keys and values {timestamp: coordinates}
            
            # Check if the previous timestamp is contained within the event
            if list(prev_event.keys())[-1] == str(prev_ts):
                
                # Check if the two events also overlapse in space, and not only in time  
                if isOverlaping(current_event, event_area[name_prev][str(prev_ts)]):                # True, if they overlap in area 
                    #event_area[name_prev][str(current_ts)] = copy.deepcopy(current_event)      #BRAUCHTS DAS???  # append the {'timestamp': coordinates}-pair of the next timestep to this event  
                    existingEvent       = True
                elif isContiguousGridPoint(current_event, event_area[name_prev][str(prev_ts)]):     # True, if it's a congiuous grid point
                    existingEvent       = True
                else:                                                                       # True, if previous timestep exists, but areas don't overlap
                    pass        
            else:   # reached, if previous timestamp is not in event. This means the event is over, and can be written to finished_events 
                finished_events.append(name_prev)
            
            if existingEvent: # enters, if current event is continuation to an already ongoing event
                updated_event_area[name_prev].update({str(current_ts): current_event})
                notBeenAttributed = False
                
        if notBeenAttributed:   # True, if current event can't be attributed to an already existing event --> write new event
            last_key = [k for k in event_area.keys()][-1]   # inefficient
            name_next = 'event' + str(int(last_key[5:len(last_key)]) + cntr)
            cntr += 1
            updated_event_area.update({name_next: {str(current_ts): current_event}})
            
            #name_next = 'event' + str(int(name_prev[5:len(name_prev)]) + cntr)
            #cntr += 1
            #updated_event_area.update({name_next: {str(current_ts): current_event}})
            
            
    finished_events_update = list(dict.fromkeys(finished_events))   # removes duplicates
    return [updated_event_area, finished_events_update]
#TODO: what if two COLs merge into one?
#TODO: get rid of daily events


#####################################################################################################
# Takes two numpy arrays, containing coordinates. Returns True, if they overlapse in area.
#
# Input:    array1                  - (n,2) numpy array
#           array2                  - (m,2) numpy array
# Output:   True or False
#####################################################################################################
def isOverlaping(array1, array2):
    for i in np.arange(len(array1)):
        for j in np.arange(len(array2)):
            if (array1[i,0] == array2[j,0]) & (array1[i,1] == array2[j,1]):     # yes, they overlap
                return True
    return False    # no, they don't overlap


#####################################################################################################
# Takes two numpy arrays, containing coordinates. Returns True, if they are contiguous gridpoints.
#
# Input:    array1                  - (n,2) numpy array
#           array2                  - (m,2) numpy array
# Output:   True or False
#####################################################################################################
def isContiguousGridPoint(array1, array2):  #TODO: do they assume 4pt-neighbourhood or 8pt???
    #for arr1_entry in array1:
    #neighbours = [] #assuming 4pt neighbourhood
    #neighbours.append([arr1_entry[0], arr1_entry[1]-1])
    #neighbours.append([arr1_entry[0], arr1_entry[1]+1])
    #neighbours.append([arr1_entry[0]-1, arr1_entry[1]])
    #neighbours.append([arr1_entry[0]+1, arr1_entry[1]])
    neighbours = getNeighbour(array1,1) #assuming 8pt neighbourhood 
    for this_neighs in neighbours:
        for this_arr2 in array2:
            if (this_neighs == this_arr2).all():
                return True

'''
#####################################################################################################
# Takes a dict of COL events and removes events that last less than two days. Input has in the format
# of the sortIntoEventArea-function output (see function documentation). Output is in the same order.
#
# Input:    event_area              - dict (output of sortIntoEventArea)
# Output:   event_area2days         - dict (same format as input)
#####################################################################################################
def deleteSingleDayEvent(event_area):
    event_area2days = {}
    event_cntr = 0
    for event_No in range(len(event_area)):     # runs trough the first layer ('eventx' layer)
        current_events = event_area['event' + str(event_No)]
        if len(current_events) > 1:   # True if it is an event lasting at least two days
            name_event = 'event' + str(event_cntr)
            event_cntr += 1     
            event_area2days.update({name_event: event_area['event'+str(event_No)]})
    return event_area2days
'''


#####################################################################################################
# Takes a dict of COL events and removes events that last less than two days. Input has in the format
# of the sortIntoEventArea-function output (see function documentation). Output is in the same format.
# Also outputs updated_finished_events, which is the same is finished_events, but without the 
# one-day-events.
#
# Input:    event_area              - dict (output of sortIntoEventArea)
#           finished_events         - list, containing str ('event n', 'event n + x', ...)
# Output:   event_area2days         - dict (same format as input)
#           return_finished_events  - list, same as finished_events, but without one-day-events
#####################################################################################################
def deleteOneDayEvents(event_area, finished_events):
    updated_finished_events = []
    for name_event in finished_events:
        if len(event_area[name_event]) == 1:    # True, if current event contains only one timestep --> One Day Event 
            if type(event_area[name_event]) is not dict:    # Checking for errors
                print('Error occured. Check name_event = ', name_event)
                sys.exit()
            del event_area[name_event]
        else:   # enters, if it is a finished, multiple-day-event
            updated_finished_events.append(name_event)  # containing names of finished, multiple-day events
    
    # event_area now only contains multiple day events, but the first-layer-keys contain gaps --> rename keys
    updated_event_area = {}
    cntr = 0
    my_keys = [key for key in event_area.keys()]    # list comprehension
    return_finished_events = []
    for event_No in range(len(event_area.keys())):
        new_name = 'event' + str(event_No)
        updated_event_area[new_name] = copy.deepcopy(event_area[my_keys[cntr]])
        if my_keys[cntr] in updated_finished_events:
            return_finished_events.append(new_name)
        cntr += 1
        
    return [updated_event_area, return_finished_events]
    
    
