#! /usr/bin/env python

"""
PCMF.py -> Plot and Compare Monitor Files

This code reads SWIM monitor file/files and plots the time-based variables (i.e. first dimension 
is 'timeDim').  Scalars are plotted as functions of time.  Profiles are plotted as functions
of rho for selected time points.  All the plots are written into a multipage pdf document.  Each
variable is plotted on a separate page.  There are summary pages of the global monitor file
attributes at the beginning of the document, and there is in index of all plots at the end of the
document.

The code takes one command-line argument -> either a  monitor file name with extension ".nc"
or a configuration file with extension ".config".  

If the file is a monitor file this code plots the scalars as functions of time (one variable 
per plot) and profiles as functions of rho for a set of time points. There are n_time_points 
(see parameter below) one curve per time point that are approximately evenly spaced in time.

If the file is a configuration file you can specify a group of monitor files to be read and
compared.  Scalars and profiles are plotted with one curve from each monitor file, one 
variable per page.  Profiles at different time points are plotted on different pages.  The
configiruation file allows you to specify other plot parameters (see discussion below).  You
can also specify a single monitor file in the config file and control the other plot
paramaters for that.

The configuration file is read by simple_config.py. It must contain: MONITOR_FiLE_NAMES and 
TIME_POINT_MODE.

It can optionally contain: PLOT_FILE_NAME, MONITOR_FILE_DIR, CURVE_LABELS 

MONITOR_FiLE_NAMES = The names of the monitor files to be plotted.  It can consist of a single
monitor file or a list of multiple monitor files.  For profile plot pages the curve labels are
different depending on whether a single file or list of files is being plotted.  If a single 
monitor file is given, curves for each of the specified time points are plotted on the page
for each variable.  If a list of monitor files is given, curves for each of the monitor files 
are plotted on the page, with the various time points on separate pages.

TIME_POINT_MODE = Specification of time points for which profiles are to be plotted.
For now there are two modes EXPLICIT or REGULAR

For TIME_POINT_MODE = REGULAR  The number of time points to be plotted, N_TIME_POINTS, must 
be specified.
 If N_TIME_POINTS = 0, only the first time point it plotted
 If N_TIME_POINTS = -1, only the last time point it plotted
 If N_TIME_POINTS <= the number of time points in the monitor file, the first, the last,
 and n-2 time points spaced as evenly as possible will be plotted

For TIME_POINT_MODE = EXPLICIT -> not implemented yet

Optional parameters are:

PLOT_FILE_NAME = Name of .pdf file where plots go.  If omitted and MONITOR_FiLE_NAMES is a
single file, the monitor file name (minus .nc extension) is used (e.g. <monitor_file>.pdf).  
If omitted and MONITOR_FiLE_NAMES is a list of monitor files, the config file name (minus .config
extension) is used (e.g. <config_file>.pdf).

MONITOR_FILE_DIR = path (unix style) to a directory containing the monitor files.  If omitted
the monitor files are assumed to be in the current working directory. 

CURVE_LABELS = A space delimited list of labels to subsitute for file names when plotting 
profiles for a list of monitor files.  When plotting profiles for a list of monitor files the 
default labels for the curves are the names of the corresponding monitor files. However this 
probably doesn't describe the differences in the simulation runs very well. So you can optionally
provide this list of substitute curve lables. There must be a one-to-one mapping between the 
labels and the monitor file names.

This script requires external modules:
    Scientific.IO.NetCDF
    plt_XY_curves - which also requires matplotlib
    simple_config

change log:
 1/5/2011
 version 1.0
 
 Version 1.1 3/30/2011 Batchelor
 Added capability to plot several monitor files with different numbers of time points.  However 
 there is a complication when plotting profiles.  Normally profiles are plotted for all monitor files
 at a fixed time.  However if the time points are different in different files how do you do this?
 This version can do this for the first and last time points (i.e. N_TIME_POINTS = 0 or -1 in the 
 config file). Hopefully these time points will line up.  If not, caveat emptor.

"""

import sys
import os

from ips_monitor.plt_XY_Curves import *
from ips_monitor.simple_config import *
#from Scientific.IO.NetCDF import *
#from scipy.io.netcdf import *
from netCDF4 import *

n_time_points = 10 # Number of time points at which to plot profiles

debug = False

#----------------------------------------------------------------------------------------------        
# Utility functions
#----------------------------------------------------------------------------------------------

def  n_evenly_spaced_integers(n, Length):

# The purpose is to generate a list of integers that can be used to index a subset of a long 
# list of length = Length, at an approximately even stride, but including the first and last
# items in the list.
# It tries to fit n evenly spaced integers in the range 0:Length-1, provided n >= Length
# If n == 0 it returns [0] i.e. it indexes the first item in the list
# If n == -1 it returns [-1] i.e. it indexes the last item in the list
# If n >= Length you are asking for more indices than would be in the list. It returns range(Length)
# If n < Length-1 returns [0, n-2 integers spaced as evenly as possible, Length-1].  Of
# course the returned list won't be exactly evenly spaced unless n divides Length.

    if (type(n) != int) or (type(Length) != int) or (Length < 1):
        print('error n_evenly_spaced_integers: arguments must be positive integers',\
              ' n = ', n, 'Length = ', Length)
              
    if n >= Length:
        return list(range(Length))
    elif n == 0:
        return [0]
    elif n == -1:
        return [-1]
    else:
        return [int((i)*float(Length-1)/(n-1)) for i in range(n)]        

#----------------------------------------------------------------------------------------------

def select_nearest_points(points, selectors):

# From list of real numbers -> 'points' find the indices in the list of the points nearest
# the values in the list -> 'selectors'.  That way you don't have to know exactly what's in
# points to get something close ot what you want.
    
    indices = []
    for target in selectors:
        diffs = [abs(x - target) for x in points]
        i = diffs.index(min(diffs))
        indices.append(i)
    return indices
    

#----------------------------------------------------------------------------------------------
# Set up
#----------------------------------------------------------------------------------------------

# get the command line -> .config or .nc file name
if len(sys.argv) != 2:
    print(' sys.argv = ', sys.argv)
    print('Usage: this script takes one command line argument -> a monitor\
           or config file name')
    raise Exception('Usage: this script takes one command line argument -> a monitor\
           or config file name')

# get the configuration -> plot file name, monitor files to compare + other stuff
if sys.argv[1][-7:] == '.config':
    config_file_name = sys.argv[1]
    config = simple_config(config_file_name)
    
    # Get the monitor file names and monitor file directory path if provided
    monitor_file_names = config.MONITOR_FILE_NAMES
    monitor_file_dir = os.getcwd()
    try:
        monitor_file_dir = config.MONITOR_FILE_DIR
    except Exception:
        print('No MONITOR_FILE_DIR provided. Using current working directory')
    
    # Determine if monitor_file_names is a single file (i.e. a string) or a list
    if type(monitor_file_names) == str:
        single_file = True
        monitor_file_names = [monitor_file_names]   # Make it into a length 1 list
    
    elif type(monitor_file_names) == list:
        single_file = False
        
        # Generate a list of lables for the curves indicating which monitor file goes with which curve.
        curve_labels = monitor_file_names       # Default is to use the file names.
        
        # If a list of curve labels is provided use that instead.
        curve_label_dict = {}
        try:
            curve_labels = config.CURVE_LABELS
            for i in range( len(monitor_file_names) ):
                curve_label_dict[ monitor_file_names[i] ] = curve_labels[i]
        except Exception:
            print('No labels provided.  Using file names')
            for i in range( len(monitor_file_names) ):
                curve_label_dict[ monitor_file_names[i] ] = monitor_file_names[i][: -3]
    
    # Get plot file name and open file for plotting
    try:
        plot_file_name = config.PLOT_FILE_NAME + '.pdf'
        given_plot_file_name = True
        print('plot_file_name = ', plot_file_name)
        if plot_file_name.strip() == '.pdf':    # PLOT_FILE_NAME parameter is present but blank
            given_plot_file_name = False
    except Exception:
        given_plot_file_name = False
    
    if given_plot_file_name == False:   # Generate default plot_file_name   
        if single_file:
            plot_file_name = monitor_file_names[0].partition('.nc')[0] + '.pdf'
        else:
            plot_file_name = config_file_name.partition('.config')[0] + '.pdf'          
        print('No plot_file_name provided.  Using default name = ', plot_file_name)

# If single monitor file, define monitor_file_names list
elif sys.argv[1][-3:] == '.nc':
    monitor_file_names = [ sys.argv[1] ]
    monitor_file_dir = os.getcwd()
    single_file = True
    plot_file_name = monitor_file_names[0][:-3] + '.pdf'
    print('plot_file_name = ', plot_file_name)
    
else:
    print('Unrecognized cammand line argument.  Should be a .nc or .config file')
    raise exception('Unrecognized cammand line argument.  Should be a .nc or .config file')

open_file_XY_Curves_Fig(plot_file_name)

all_monitor_files_dict = {}
all_monitor_var_names = []  # List of variables in all files

scalar_names = []
profile_names = []
dimensions_dict = {}
index = []
plot_count = 0


for file in monitor_file_names:
    print('Processing monitor file ', file)
    full_path = os.path.join(monitor_file_dir, file)
    monitor = Dataset(full_path, 'r', format = 'NETCDF3_CLASSIC')
    
    plot_count = plot_count +1
    monitor_dim_names = list(monitor.dimensions.keys())
    monitor_dims = monitor.dimensions
    monitor_var_names = list(monitor.variables.keys())

    # Compile list of all variables in all files. Add new variable first time it appears.
    # Determine if scalar or profile.  Record variable dimensions.
    for var_name in monitor_var_names:
        if var_name not in all_monitor_var_names: 
            all_monitor_var_names.append(var_name)
            if var_name == 'time':   # skip time variable
                continue
            var_shape = monitor.variables[var_name].shape
            var_dimensions = monitor.variables[var_name].dimensions

            if debug:
                print(var_name, ' shape =  ', var_shape)
                print(var_name, ' dim names =  ', monitor.variables[var_name].dimensions)
                    
            # Only consider time based variables i.e. first dim == timeDim
            first_dim = var_dimensions[0]
            if first_dim == 'timeDim':
                if len(var_shape) == 1:     # i.e. is scalar
                    scalar_names.append(var_name)
                elif len(var_shape) == 2:   # i.e. is profile
                    profile_names.append(var_name)

            dimensions_dict[var_name] = getattr(monitor.variables[var_name], 'units')
            
    # Global attributes defined in this monitor file
    Global_label = monitor.Global_label
    RunID = monitor.RunID
    tokamak_id = monitor.tokamak_id
    shot_number = monitor.shot_number
    
    global_attributes = [ ['Monitor file = ', file],\
                          ['Global label = ', Global_label], ['Run ID = ', RunID],\
                          ['tokamak id = ', tokamak_id], ['shot number = ', str(shot_number)] ]
    
    if debug:
        print('Global_label = ', Global_label)
        print('RunID = ', RunID)
        print('tokamak_id = ', tokamak_id)
        print('shot_number = ', shot_number)
        print(' ')
        monitor_dim_names.sort()
        print('monitor_dims = ', monitor_dim_names)
        print(' ')
        monitor_var_names.sort()
        print('monitor_var_names = ', monitor_var_names)
        print(' ')
        #monitor_dims.sort()
        print('monitor_dims = ', monitor_dims)
        print(' ')
    
    summary ={'file_name':file , 'global_attributes': global_attributes}
    plot_summary(summary)

    index.append([plot_count, file])
        
    # Get time points for this monitor file
    time = monitor.variables['time'][:]
    print('time = ', time)

    this_monitor_file_dict = {'globals': global_attributes, 'time': time,\
                              'variables': monitor_var_names, 'file': monitor}
    all_monitor_files_dict[file] = this_monitor_file_dict

#----------------------------------------------------------------------------------------------
# Plot all scalars versus time, if there is more than one time point
#----------------------------------------------------------------------------------------------

if len(time) > 1:   
    print('Plotting scalars \n')
    
    for var_name in scalar_names:
        if debug:
            print('var_name = ', var_name)
    
        title = var_name
        xlabel = 'time (sec)'
        
        ylabel = var_name
        units = dimensions_dict[var_name]
        if units != ' ' and units != '':
            ylabel = var_name + ' (' + units + ') '
    
        # Add a curve from each monitor file that contains this variable
        curve_list = []
        for file_name in monitor_file_names: 
            this_monitor_file_dict = all_monitor_files_dict[file_name]
            monitor = this_monitor_file_dict['file']
            time = this_monitor_file_dict['time']
            variables = this_monitor_file_dict['variables']
            
            # Check if var_name is in this monitor file.
            if var_name in variables:
                var = monitor.variables[var_name]
                y = var[:]
                
                if debug:
                    print(ylabel)
                    print('type(', var_name, ') = ', type(y))
                    print(var, ' = ', y)
                    
                if single_file:
                    curve_list.append( XY_curve(time, y) )
                else:
                    curve_list.append( XY_curve(time, y, label = curve_label_dict[file_name]) )
    
        plot1 = XY_Curves_Fig(curve_list, title, xlabel, ylabel)
        plot_XY_Curves_Fig(plot1)
        plot_count = plot_count + 1
        index.append([plot_count, title])
    
#----------------------------------------------------------------------------------------------
# Plot all profiles versus rho, at an evenly spaced set of time points -> n_time_points
#----------------------------------------------------------------------------------------------
print('Plotting profiles \n')

# get the configuration -> time points
if sys.argv[1][-7:] == '.config':
    
    # Get config parameters defining time points to be plotted
    time_point_mode = config.TIME_POINT_MODE
    
    if time_point_mode == 'REGULAR':
        n_time_points = int(config.N_TIME_POINTS)
        time_point_indices =  n_evenly_spaced_integers(n_time_points, len(time))
    
    if time_point_mode == 'EXPLICIT':
        requested_time_points_str = config.TIME_POINTS
        requested_time_points = list(map(float, requested_time_points_str))
        time_point_indices = select_nearest_points(time, requested_time_points)

# If single monitor file, use default n_time_points defined at top
elif sys.argv[1][-3:] == '.nc':
        time_point_indices = n_evenly_spaced_integers(n_time_points, len(time))

print('n_time_points = ', n_time_points, '   time_point_indices = ', time_point_indices)
print('time_points = ', [time[i] for i in time_point_indices])

for var_name in profile_names:
    if debug:
        print('var_name = ', var_name)

    xlabel = 'rho'    
    ylabel = var_name
    units = dimensions_dict[var_name]
    if units != ' ' and units != '':
        ylabel = var_name + ' (' + units + ') '

    # skip rlim and zlim since they aren't profiles
    if var_name == 'rlim' or var_name == 'zlim':
        continue
    if debug:
        print('var_name = ', var_name)

    curve_list = []

    # Add a plot of this variable for each time point
    for i in time_point_indices:

        #Generate title
        title = var_name
        if not single_file:
            time_label = 't = ' + str(time[i])
            title = var_name + '    ' + time_label
        

        # Add a curve from each monitor file that contains this variable
        if not single_file:
            curve_list = []
        for file_name in monitor_file_names: 
            this_monitor_file_dict = all_monitor_files_dict[file_name]
            monitor = this_monitor_file_dict['file']
            variables = this_monitor_file_dict['variables']
            
            # Check if var_name is in this monitor file.
            if var_name in variables:
                var = monitor.variables[var_name]
                y = var[:]
                
                if debug:
                    print(ylabel)
                    print('type(', var, ') = ', type(y))
                    print(var, ' = ', y)
                
                # Get the grid variable from the dimension name    
                # Also determine if this is a zone based grid, indicated by prefix "dm1_n', if so 
                # delete the "dm1_n" and append "_zone" so we can find it
                grid_name = var.dimensions[1]
                if grid_name[0:5] == 'dm1_n':
                    grid_name = grid_name[5:] + '_zone'
                    grid = monitor.variables[grid_name][:]           
                elif grid_name[0:5] == 'dim_n':
                    grid_name = grid_name[5:]
                    grid = monitor.variables[grid_name][:]
                else:
                    message = 'This dimension name -> ' + grid_name + ' does not correspond to a grid'
                    print(message)
                    raise Exception(message)
            
                if debug:
                    print('grid = ', grid)
                    
                if single_file:
                    lbl = 't = ' + str(time[i])
                    curve_list.append( XY_curve(grid, y[i], label = lbl))                
                else: 
                    curve_list.append( XY_curve(grid, y[i], label = curve_label_dict[file_name]) )  
    
        if not single_file:     # Generate the plot just outside the file_name loop
            plot1 = XY_Curves_Fig(curve_list, title, xlabel, ylabel)
            plot_XY_Curves_Fig(plot1)
            plot_count = plot_count + 1
            index.append([plot_count, title])

    if single_file:     # Generate the plot just outside the time_point_indices loop
        plot1 = XY_Curves_Fig(curve_list, title, xlabel, ylabel)
        plot_XY_Curves_Fig(plot1)
        plot_count = plot_count + 1
        index.append([plot_count, title])

#----------------------------------------------------------------------------------------------
# Finalize
#----------------------------------------------------------------------------------------------

plot_index(index, 1)
    
close_file_XY_Curves_Fig()

if debug:
    print('index = ', index)
