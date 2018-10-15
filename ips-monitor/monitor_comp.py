#! /usr/bin/env python

"""
monitor_comp.py 10-11-2018

Version 5.4 runs PCMF to generate pdf plot of plasma state and pushes it out to the w3 
directory.

Version 5.3 picks up the simulation config file and pushes it out to the w3 directory.
This is done in the init

Version 5.2 adds a new variable V_loop(rho) and removes requirement for having an ElViz
template file.  Nobody seems to be using ElViz anymore, but until now the component
looked for a template .xml file (specified in the simulation config file) and would 
throw an exception if it was missing.  No more.  Users can still use ElViz if they like
but they must access it directly.  Users can, and probably should eliminate the
template file from the [MONITOR] section of the simulation config file.  Many config files
now have a line "INPUT_FILES = basic_time_traces.xml", The framework would throw an
exception if the .xml file weren't there.

Version 5.1 adds three new veriables Z_axis, P_nuc_FUS and Q_nuc_FUS.  This takes the 
instantaneous fusion power production, P_nuc_FUS, to be 5 times the instantaneous plasma 
heating by fusion products.  It might not be such a good approximation for a machine 
smaller than ITER. The calculation of P_nuc_FUS is sort of interesting because it calls 
calculate_MonitorVariable( ) recursively.  And Q_nuc_FUS is interesting because the 
monitorDefinition dictionary gets modified in the init_monitor_file() function.  See notes 
from 4/2/12 below.

Version 5 is a major change from previous versions in that in addition to keeping power
and current profiles as zone based quantities (Watts/zone or Amps/zone) it also gives
power and current densities (Watts/m^3 and Amps/m^2).  Plus it adds cumulative 
profiles for these quantities (Watts < rho and Amps < rho).  It uses the netCDF4.py
module as did the last version, monitor_comp_4. The name convention is for example:
    Pe_NB_dens(rho) = NB electron heating power density profile in W/m^3
    Pe_NB(rho) = NB electron heating power deposition profile in W/zone
    Pe_icrf_cum(rho) = NB electron heating power deposited inside rho in W
    J_NB(rho) = Beam driven current density profile in A/m^2    
    I_NB(rho) = Beam driven current profile in A/zone   
    I_NB_cum(rho) = Beam driven current inside rho in A

The user selects which of these sets of quantities (densities, zone based, cumulative) to
include in the monitor file by setting the following variables in the [MONITOR] section 
of the simulation config file:
    INCLUDE_DENSITY_BASED_PROFILES = True/False
    INCLUDE_ZONE_BASED_PROFILES = True/False
    INCLUDE_CUMULATIVE_PROFILES = True/False
The densities are in some sense more useful in that profiles defined on different
radial grids can be quantitatively compared, whereas profiles based on zones cannot.
Therefore the default is to include density based profiles and cumulative profiles,
but this behavior can be changed by setting monitor component configuration parameters
INCLUDE_DENSITY_BASED_PROFILES, INCLUDE_ZONE_BASED_PROFILES or
INCLUDE_CLUMULATIVE_PROFILES to True/Flase in the simulation config file.
    
    
A Python component to replace the previous IPS monitor component written in f90 (i.e.
monitor_comp_4.f90).

The intention is that this component be partially configurable from the SWIM .congfig
file. That is one can select from a list of quantities to monitor.  N.B. This 
.config file configuration has now been partially implemented in the guise of
selecting between density, zone based, and cumulative profiles.  So far there has been
no demand for selecting individual variables in the config file.

This component has
to know how to calculate the things it's requested to monitor, so there is at least some  
custom programing required for each quantity that can be monitored.  To add a new
variable to monitor there are three customizations that must be done:

1) The name of the new variable, in the form you want the lable displayed, must be
   added to the 'requestedVars' list below.

2) There must be an entry in the 'monitorDefinition' dictionary giving the name and a 
   list of properties for the monitor variable.  See discussion below on
   'monitorDefinition'.

3) There must be an 'if' block in the calculate_monitorVariable function that calculates
   the value of the monitor variable from the data it depends upon.

This script exports 4 SWIM component functions:
init - defines the netcdf monitor file for monitor data then calls step to insert the
       initial data (t = t0).  Also saves some internal data needed to restart the
       monitor component (monitorVars and ps_VarsList).  These are are pickled to file
       'monitor_restart'
restart - unpickles the internal monitor data previously saved in file 'monitor_restart'
          and loads these global variables.
step - writes monitor data for current time step to the monitor netcdf file and saves it
finalize - sets the global netcdf attribute to tell ELViz to stop watching this
           monitor file

change log:
 1/21/09
 version 4.1 adds the variable kind '2D' and implements one such variable 'Psi(RZ)'
 2/20/09
 version 4.2 adds variables for EC and NB components plus a few other state variables
 3/17/09
 version 4.3 adds variables for Fusion Fast Ion species
 4/7/09
 version 4.4 adds <Te>, <Ti>, <ne>, I_BS, Pe_OH, q(95), and li(3)
 8/19/09
 version 4.5 adds global attributes Global_label, RunID, tokamak_id, shot_number
             adds variables I_OH, Zeff, triang, elong, rlim, and zlim
             To accommodate  rlim, and zlim also added a new kind -> scalar list
 3/22/10
 version 4.6 adds variables for Lower Hybrid heating and current drive
  
 4/11/10
 version 4.7 adds restart function and saves internal state to file monitor_restart

 5/24/10
 version 4.8 adds checkpoint function

 10/7/10
 version 4.9 adds new variables N_GW, beta_th, beta_N, n_H, n_D, n_T, n_He4 
             Also sets a limit (presently = 10) on the maximum vaule of q(0), q(95) and 
             q(rho).  This is so the plots won't get so squeezed

1/30/11 4.10 Replaces the obsolete Scientific.IO.NetCDF.py module with netCDF4.py.  The 
             netCDF4.py module is more complete and much better supported than the Scientific
             version. It still uses the netcdf3 interface but that is also supported by the 
             netCDF4.py module.  To use this you need netCDF4 in your python path.
        
3/13/12 5.0  Adds power and current density profiles and cumulative profiles in addition
             to the zone based profiles of previous versions

4/2/12 5.1  Adds nuclear fusion power production P_nuc_FUS and Q_nuc_FUS where
             Q_nuc_FUS = P_nuc_FUS/Total heating power.  Since various heating sources
             may or may not be used in a given run and therefore may not be in the plasma 
             state, it is a little complicated to specify the plasma state dependencies
             for Q_nuc_FUS.  See working note below.  Also note that the calculation of
             P_nuc_FUS is simplistic for now.  It's just taken to be 17.6/3.5 ~5 times the
             electron and ion heating by fusion alpha's.  This is ok in steady state and
             for a machine with good alpha confinement, but probably not good for a 
             machine smaller than ITER.
             
             Also Z_axis is added.  R_axis has been there from the beginning.

11/5/2013 5.2  Adds new variable V_loop(rho).  Eliminates reference to ElViz template file
              i.e. no longer requires template file in inputs, so doesn't push a template
              file to the w3 directory.    Eliminated test for monitor component
              version number in monitor definition metadata. (So as of now there is no 
              metadata.)  It made some sense when I thought monitor variables might be 
              specified from the config file. But no more.

1/24/2014	Changed the test in function area_grid_interp() so that the validity test for
			grids has a tolerance of 1.e-12.  The aorsa component was crashing the 
			component because the the upper value of rho_icrf was not exactly 1.0.  So 
			now the first and last rho values have to be within this tolerance of 0.0
			and 1.0 respectively.  I only hope there are not other tests for exactly 0.0
			and 1.0.

2/24/2014   Version 5.3 picks up the simulation config file and pushes it out to the w3 
            directory.  The portal run_id is prepended to the file name, like with the
            monitor file.  Also cleaned up a few bogus references to monitor_comp_4 in
            print statements.
            
1/4/2018    Commented out reference to PORTAL run_id due to demise of SWIM PORTAL

1/30/2018   Replacing defunct PORTAL run_id with datetime() to distinguish between runs
"""

# Note (4/2/12)
# Some monitor variables, e.g. Q_nuc_FUS, make sense with or without some of the
# plasma state variables present.  For example Q_nuc_FUS is well defined even if some
# of the heating sources (e.g. ECH) are not present.  So need to modify the dependency list
# in the monitorDefinition dictionary for these variables based on what is actually
# in the plasma state.  In the initial specification of the monitorDefinition dictionary
# the only heating dependency listed is Ohmic.  But in the init_monitor_file() function,
# after the initial plasma state is read, there is coding to look at the plasma state and
# add as dependencies power_ec, power_lh, power_ic, and power_nbi if they are present.


# Note (10/31/10)  It is a little complicated to access some items in the plasma state
# because their names appear in plasma state lists rather than at the top level of the 
# plasma state. The prime example is thermal ion species whose names appear in the plasma
# state list 'S_name'.  For example monitor variable n_H(rho) depends on ps varible 'ns'.  We
# have to look in 'S_name' to see if an H species is actually defined, and if it is we have to
# find what index in 'ns' corresponds to it.  The way I'm handling this is to add yet another 
# dictionary containing the monitor variable name as a key and pointing to a list 
# containing the plasma state labeled list name and and the item name(s)(the variable could
# depend on more than one item in the list.  The new dictionary is called 
# 'labeled_ps_quantities'. For example for n_H(rho) the relevant entry in 
# labeled_ps_quantities dictionary is:
# 'n_H(rho)':['S_name', 'H ']
# To find out if there is a thermal H species we must look in ps list 'S_name' and see if
# label 'H' is present.  A utility function is provided to do this check and return the index
# of 'H' in the list.

import sys
import os
import subprocess
import shutil
import pickle
import time
import datetime

from  component import Component
from get_IPS_config_parameters import get_component_param

# Import the necessary Numeric and netCDF modules
from netCDF4 import *
from numpy import *


# ------------------------------------------------------------------------------
#
# Global definitions
#
# ------------------------------------------------------------------------------


monitor_comp_version = '5.3'
debug = False
include_density_based_profiles = True
include_zone_based_profiles = False
include_cumulative_profiles = True

monitor_fileName = 'monitor_file.nc'
pdf_fileName = 'monitor_file.pdf'

# List of requested variables to monitor (if dependencies are satisfied)
# The list just below is the default containing everything.  In the component it can
# be overwritten with the configuration file data (some day maybe).

requestedVars = [
    # Thermal species data
        'Te(0)',
        'Te(rho)',
        'Te_ave',
        'Ti(0)',
        'Ti_ave',
        'Ti(rho)',
        'ne(0)',
        'ne_ave',
        'ne(rho)',
        'n_H(0)',
        'n_H_ave',
        'n_H(rho)',
        'n_D(0)',
        'n_D_ave',
        'n_D(rho)',
        'n_T(0)',
        'n_T_ave',
        'n_T(rho)',
        'n_He4(0)',
        'n_He4_ave',
        'n_He4(rho)',
        'I_BS_total',
        'Pe_OH_total',
        'I_OH_total',
        'Zeff(rho)',
        'Zeff(0)',
        'Zeff_ave',
    # equilibrium data
        'P_eq(rho)',
        'q_eq(rho)',
        'q(0)',
        'q(95)',
        'li(3)',
        'R_axis',
        'Z_axis',
        'triang(95)',
        'elong(95)',
        'Psi(RZ)',
        'I_plasma',
        'rlim',
        'zlim',
        'N_GW',
        'beta_th',
        'beta_N',
    # other state data
        'Vsurf',
        'V_loop(rho)',
    # ech power deposition
        'power_EC',
        'Pe_ecrf_total',
        'I_ecrf_total',
        # Lower Hybrid power deposition
        'power_LH',
        'Pe_LH_total',
        'Pi_LH_total',
        'I_LH_total',
    # fusion product fast ion
        'ni_FUS(rho)',
        'Pe_FUS_total',
        'Pi_FUS_total',
        'I_FUS_total',
        'Eperp_FUSI(rho)',
        'Epll_FUSI(rho)',
        'P_nuc_FUS',
        'Q_nuc_FUS',
     # icrf power deposition
        'power_IC',
        'Pe_icrf_total',
        'Pi_icrf_total',
        'I_icrf_total',
     # minority species data
        'nmin_icrf(rho)',
        'Pmin_e_total',
        'Pmin_i_total',
        'Eperp_mini(rho)',
        'Epll_mini(rho)',
      # neutral beam power deposition
        'power_NB',
        'ni_NB(rho)',
        'Pe_NB_total',
        'Pi_NB_total',
        'I_NB_total',
        'Eperp_NBI(rho)',
        'Epll_NBI(rho)'
]

density_based_profiles = [
    # Thermal species data
        'J_BS(rho)',
        'Pe_OH_dens(rho)',
        'J_OH(rho)',
    # equilibrium data
        'J_plasma(rho)',
    # ech power deposition
        'Pe_ecrf_dens(rho)',
        'J_ecrf(rho)',
    # Lower Hybrid power deposition
        'Pe_LH_dens(rho)',
        'Pi_LH_dens(rho)',
        'J_LH(rho)',
    # fusion product fast ion
        'Pe_FUS_dens(rho)',
        'Pi_FUS_dens(rho)',
        'J_FUS(rho)',
     # icrf power deposition
        'Pe_icrf_dens(rho)',
        'Pi_icrf_dens(rho)',
        'J_icrf(rho)',
     # minority species data
        'Pmin_e_dens(rho)',
        'Pmin_i_dens(rho)',
      # neutral beam power deposition
        'Pe_NB_dens(rho)',
        'Pi_NB_dens(rho)',
        'J_NB(rho)'
    ]

zone_based_profiles = [
    # Thermal species data
        'I_BS(rho)',
        'Pe_OH(rho)',
        'I_OH(rho)',
    # equilibrium data
        'I_plasma(rho)',
    # ech power deposition
        'Pe_ecrf(rho)',
        'I_ecrf(rho)',
        # Lower Hybrid power deposition
        'Pe_LH(rho)',
        'Pi_LH(rho)',
        'I_LH(rho)',
    # fusion product fast ion
        'Pe_FUS(rho)',
        'Pi_FUS(rho)',
        'I_FUS(rho)',
     # icrf power deposition
        'Pe_icrf(rho)',
        'Pi_icrf(rho)',
        'I_icrf(rho)',
     # minority species data
        'Pmin_e(rho)',
        'Pmin_i(rho)',
      # neutral beam power deposition
        'Pe_NB(rho)',
        'Pi_NB(rho)',
        'I_NB(rho)'
    ]

cumulative_profiles = [
        'I_BS_cum(rho)',
        'Pe_OH_cum(rho)',
        'I_plasma_cum(rho)',
        'I_OH_cum(rho)',
        'Pe_ecrf_cum(rho)',
        'Pe_LH_cum(rho)',
        'Pi_LH_cum(rho)',
        'I_LH_cum(rho)',
        'Pe_FUS_cum(rho)',
        'Pi_FUS_cum(rho)',
        'I_FUS_cum(rho)',
        'Pe_icrf_cum(rho)',
        'Pi_icrf_cum(rho)',
        'I_icrf_cum(rho)',
        'Pe_NB_cum(rho)',
        'Pi_NB_cum(rho)',
        'I_NB_cum(rho)'
#        'Area(rho)',
#        'Volume(rho)'
        ]

# List of variables to monitor (which dependencies are satisfied)
monitorVars = []

# List of grids needed to plot monitor variables
monitorGrids = []

# List of needed Plasma State variables needed to calculate monitor variables
ps_VarsList = []

# List of non-grid dimensions needed for other variables - e.g. scalar lists
monitorDims = []

# Dictionary of Plasma State dependencies for each variable to be monitored:
# For the time being assume that all dependencies are in the Plasma State and that any
# dimensions of the monitor variables will concide with a dimension in the Plasma State 
# (other than the unlimited time dimension). The dictionary entries have the form:
# 'variable name' : [ list of properties ].
# The list of properties is of length 3 for scalars and 4 for profiles, lists, or 2D. 
# it contains:
# (1) 'kind
# (2) 'units'
# (3) [list of ps variables upon which the monitor variable depends]
# (4) [list of grids for plotting or dimensions for lists  (note: this list is absent for 
#      scalars)]
#
# kind = one of 'S' -> scalar, 'P' -> profile, '2D' -> {X,Y} array (actually usually {R,Z}
# or 'SL' -> list of scalars.
# If kind = 'P' there must be at least one grid specified for plotting
#
# if the name of a plotting grid has '~' prepended it's assumed that the variable is zone-
# centered in this grid and a new grid is generated, <grid>_zone with size = size(<grid> - 1
# If kind = 'SL' there must be one dimension specified as the index of the list
#
# Note about order of grids in 2D variables:  In plasma state things called something like
# Psi(R,Z) are actually dimensioned PsiRZ(dim_nz, dim_nr).  I don't know if this applies to
# all 2D things in the plasma state, but it's advisable to check this in plasma state when 
# setting up 2D variables in the monitor definition section and get the right order.  If
# you don't you will get 'shapes not aligned errors when trying to set the netcdf monitor 
# variables to the calculated values.
#
# Note about dimensions in the Plasma State netcdf file:  The ps netcdf file does not 
# explicitly contain variables associated with the dimensions in the plasma state. Rather these
#  appear as netcdf dimensions.  For example there is no nc variable 'nrho' but there is an nc
# dimension 'dim_nrho'.  For safety, and simplicity I have not relied on the dimension values
# but work from the lengths of the variables given by 'len(variable)'


monitorDefinition = {           # metadata
            'monitor_comp_metaData': [],
                                # Thermal species data
            'Te(0)':['S', 'keV', ['Ts'] ],
            'Te_ave':['S', 'keV', ['Ts', 'rho', 'rho_eq', 'vol'] ],
            'Ti(0)':['S', 'keV', ['Ti'] ],
            'Ti_ave':['S', 'keV', ['Ti', 'rho', 'rho_eq', 'vol'] ],
            'Te(rho)':['P', 'keV', ['Ts'],['~rho'] ],
            'Ti(rho)':['P', 'keV', ['Ti'],['~rho'] ],
            'ne(0)':['S', 'm^-3', ['ns'] ],
            'ne_ave':['S', 'm^-3', ['ns', 'rho', 'rho_eq', 'vol'] ],
            'ne(rho)':['P', 'm^-3', ['ns'],['~rho'] ],
            'n_H(0)':['S', 'm^-3', ['ns', 'S_name'] ],
            'n_H_ave':['S', 'm^-3', ['ns', 'S_name', 'rho', 'rho_eq', 'vol'] ],
            'n_H(rho)':['P', 'm^-3', ['ns', 'S_name'],['~rho'] ],
            'n_D(0)':['S', 'm^-3', ['ns', 'S_name'] ],
            'n_D_ave':['S', 'm^-3', ['ns', 'S_name', 'rho', 'rho_eq', 'vol'] ],
            'n_D(rho)':['P', 'm^-3', ['ns', 'S_name'],['~rho'] ],
            'n_T(0)':['S', 'm^-3', ['ns', 'S_name'] ],
            'n_T_ave':['S', 'm^-3', ['ns', 'S_name', 'rho', 'rho_eq', 'vol'] ],
            'n_T(rho)':['P', 'm^-3', ['ns', 'S_name'],['~rho'] ],
            'n_He4(0)':['S', 'm^-3', ['ns', 'S_name'] ],
            'n_He4_ave':['S', 'm^-3', ['ns', 'S_name', 'rho', 'rho_eq', 'vol'] ],
            'n_He4(rho)':['P', 'm^-3', ['ns', 'S_name'],['~rho'] ],
            'J_BS(rho)':['P', 'A m^-2', ['curr_bootstrap','rho', 'rho_eq', 'area'],['rho_eq'] ],
            'I_BS(rho)':['P', 'A', ['curr_bootstrap'],['~rho'] ],
            'I_BS_cum(rho)':['P', 'A', ['curr_bootstrap'],['~rho'] ],
            'I_BS_total':['S', 'A', ['curr_bootstrap'] ],
            'J_OH(rho)':['P', 'A m^-2', ['curr_ohmic','rho', 'rho_eq', 'area'],['rho'] ],
            'I_OH(rho)':['P', 'A', ['curr_ohmic'],['~rho'] ],
            'I_OH_cum(rho)':['P', 'A', ['curr_ohmic'],['~rho'] ],
            'I_OH_total':['S', 'A', ['curr_ohmic'] ],
            'Pe_OH_dens(rho)':['P', 'W m^-3', ['pohme','rho', 'rho_eq', 'vol'],['rho'] ],
            'Pe_OH(rho)':['P', 'W', ['pohme'],['~rho'] ],
            'Pe_OH_cum(rho)':['P', 'W', ['pohme'],['~rho'] ],
            'Pe_OH_total':['S', 'W', ['pohme'] ] ,
            'Zeff(rho)':['P', ' ', ['Zeff'],['~rho'] ],
            'Zeff(0)':['S', ' ', ['Zeff'] ],
            'Zeff_ave':['S', ' ', ['Zeff', 'rho', 'rho_eq', 'vol'] ],

                                # equilibrium data
            'P_eq(rho)':['P', 'Pa',['P_eq'],['rho_eq'] ],
            'q_eq(rho)':['P', ' ',['q_eq'],['rho_eq'] ],
            'q(0)':['S', ' ', ['q_eq'] ],
            'q(95)':['S', ' ',['q_eq', 'psipol','rho_eq'] ],
            'li(3)':['S', ' ',['curt', 'q_eq', 'grho2r2i', 'phit', 'R_axis', 'vol'] ],
            'R_axis':['S', 'm', ['R_axis'] ],
            'Z_axis':['S', 'm', ['Z_axis'] ],
            'triang(95)':['S', ' ',['triang', 'psipol','rho_eq', 'rho_eq_geo'] ],
            'elong(95)':['S', ' ',['elong', 'psipol','rho_eq', 'rho_eq_geo'] ],
            'Psi(RZ)':['2D', 'Wb/rad', ['PsiRZ'], ['Z_grid', 'R_grid'] ],
            'rlim':['SL', 'm', ['rlim'], ['dim_num_rzlim'] ],
            'zlim':['SL', 'm', ['zlim'], ['dim_num_rzlim'] ],
            'J_plasma(rho)':['P', 'A m^-2', ['curt', 'rho_eq', 'area'],['rho_eq'] ],
            'I_plasma(rho)':['P', 'A', ['curt'],['~rho_eq'] ],
            'I_plasma_cum(rho)':['P', 'A', ['curt'],['rho_eq'] ],
            'I_plasma':['S', 'A', ['curt'] ],
            'N_GW':['S', ' ', ['ns', 'rho', 'rho_eq_geo', 'R_surfMin', 'R_surfMax','curt'] ],
            'beta_th':['S', ' ',['P_eq', 'rho', 'rho_eq','vol', 'B_axis_vac'] ],
            'beta_N':['S', ' ', ['P_eq', 'rho', 'rho_eq','vol', 'B_axis_vac',\
                                 'R_surfMax', 'R_surfMin', 'curt'] ],

                                # other state data
            'Vsurf':['S', 'Volts', ['vsur'] ],
            'V_loop(rho)':['P', 'Volts', ['V_loop'],['rho'] ],
                                # ech power deposition
            'power_EC':['S', 'W', ['power_ec'] ],
            'Pe_ecrf_dens(rho)':['P', 'W m^-3', ['peech','rho_ecrf', 'rho_eq', 'vol'],['rho_ecrf'] ],
            'Pe_ecrf(rho)':['P', 'W', ['peech'],['~rho_ecrf'] ],
            'Pe_ecrf_cum(rho)':['P', 'W', ['peech'],['~rho_ecrf'] ],
            'Pe_ecrf_total':['S', 'W', ['peech'] ] ,
            'J_ecrf(rho)':['P', 'A m^-2', ['curech','rho_ecrf', 'rho_eq', 'area'],['rho_ecrf'] ],
            'I_ecrf(rho)':['P', 'A', ['curech'],['~rho_ecrf'] ],
            'I_ecrf_cum(rho)':['P', 'A', ['curech'],['~rho_ecrf'] ],
            'I_ecrf_total':['S', 'A', ['curech'] ],
                                # LH power deposition
            'power_LH':['S', 'W', ['power_lh'] ],
            'Pe_LH_dens(rho)':['P', 'W m^-3', ['pelh','rho_lhrf', 'rho_eq', 'vol'],['rho_lhrf'] ],
            'Pe_LH(rho)':['P', 'W', ['pelh'],['~rho_lhrf'] ],
            'Pe_LH_cum(rho)':['P', 'W', ['pelh'],['~rho_lhrf'] ],
            'Pe_LH_total':['S', 'W', ['pelh'] ] ,
            'Pi_LH_dens(rho)':['P', 'W m^-3', ['pilh','rho_lhrf', 'rho_eq', 'vol'],['rho_lhrf'] ],
            'Pi_LH(rho)':['P', 'W', ['pilh'],['~rho_lhrf'] ],
            'Pi_LH_cum(rho)':['P', 'W', ['pilh'],['~rho_lhrf'] ],
            'Pi_LH_total':['S', 'W', ['pilh'] ] ,
            'J_LH(rho)':['P', 'A m^-2', ['curlh','rho_lhrf', 'rho_eq', 'area'],['rho_lhrf'] ],
            'I_LH(rho)':['P', 'A', ['curlh'],['~rho_lhrf'] ],
            'I_LH_cum(rho)':['P', 'A', ['curlh'],['~rho_lhrf'] ],
            'I_LH_total':['S', 'A', ['curlh'] ],
                                # fusion product fast ion
            'ni_FUS(rho)':['P', 'm^-3', ['nfusi'],['~rho_fus'] ],
            'Pe_FUS_dens(rho)':['P', 'W m^-3', ['pfuse','rho_fus', 'rho_eq', 'vol'],['rho_fus'] ],
            'Pi_FUS_dens(rho)':['P', 'W m^-3', ['pfusi','rho_fus', 'rho_eq', 'vol'],['rho_fus'] ],
            'Pe_FUS(rho)':['P', 'W', ['pfuse'],['~rho_fus'] ],
            'Pi_FUS(rho)':['P', 'W', ['pfusi'],['~rho_fus'] ],
            'Pe_FUS_cum(rho)':['P', 'W', ['pfuse'],['~rho_fus'] ],
            'Pi_FUS_cum(rho)':['P', 'W', ['pfusi'],['~rho_fus'] ],
            'Pe_FUS_total':['S', 'W', ['pfuse'] ],
            'Pi_FUS_total':['S', 'W', ['pfusi'] ],
            'J_FUS(rho)':['P', 'A m^-2', ['curfusn','rho_fus', 'rho_eq', 'area'],['rho_fus'] ],
            'I_FUS(rho)':['P', 'A', ['curfusn'],['~rho_fus'] ],
            'I_FUS_cum(rho)':['P', 'A', ['curfusn'],['~rho_fus'] ],
            'I_FUS_total':['S', 'A', ['curfusn'] ],
            'Eperp_FUSI(rho)':['P', 'keV', ['eperp_fusi'],['~rho_fus'] ],
            'Epll_FUSI(rho)':['P', 'keV', ['epll_fusi'],['~rho_fus'] ],
            'P_nuc_FUS':['S', 'W', ['pfuse','pfusi','rho_fus', 'rho_eq', 'vol'] ],
            'Q_nuc_FUS':['S', '', ['pfuse','pfusi','rho_fus', 'rho_eq', 'vol', 'pohme'] ],
                                # icrf power deposition
            'power_IC':['S', 'W', ['power_ic'] ],
            'Pe_icrf_dens(rho)':['P', 'W m^-3', ['picrf_totals','rho_icrf', 'rho_eq', 'vol'],['rho_icrf'] ],
            'Pi_icrf_dens(rho)':['P', 'W m^-3', ['picth','rho_icrf', 'rho_eq', 'vol'],['rho_icrf'] ],
            'Pe_icrf(rho)':['P', 'W', ['picrf_totals'],['~rho_icrf'] ],
            'Pi_icrf(rho)':['P', 'W', ['picth'],['~rho_icrf'] ],
            'Pe_icrf_cum(rho)':['P', 'W', ['picrf_totals'],['~rho_icrf'] ],
            'Pi_icrf_cum(rho)':['P', 'W', ['picth'],['~rho_icrf'] ],
            'Pe_icrf_total':['S', 'W', ['picrf_totals'] ],
            'Pi_icrf_total':['S', 'W', ['picth'] ],
            'J_icrf(rho)':['P', 'A m^-2', ['curich','rho_icrf', 'rho_eq', 'area'],['rho_icrf'] ],
            'I_icrf(rho)':['P', 'A', ['curich'],['~rho_icrf'] ],
            'I_icrf_cum(rho)':['P', 'A', ['curich'],['~rho_icrf'] ],
            'I_icrf_total':['S', 'A', ['curich'] ],
                                # minority species data
            'nmin_icrf(rho)':['P', 'm^-3', ['nmini'],['~rho_icrf'] ],
            'Pmin_e_dens(rho)':['P', 'W m^-3', ['pmine','rho_icrf', 'rho_eq', 'vol'],['rho_icrf'] ],
            'Pmin_i_dens(rho)':['P', 'W m^-3', ['pmini','rho_icrf', 'rho_eq', 'vol'],['rho_icrf'] ],
            'Pmin_e(rho)':['P', 'W', ['pmine'],['~rho_icrf'] ],
            'Pmin_i(rho)':['P', 'W', ['pmini'],['~rho_icrf'] ],
            'Pmin_e_total':['S', 'W', ['pmine'] ],
            'Pmin_i_total':['S', 'W', ['pmini'] ],
            'Eperp_mini(rho)':['P', 'keV', ['eperp_mini'],['~rho_icrf'] ],
            'Epll_mini(rho)':['P', 'keV', ['epll_mini'],['~rho_icrf'] ],
                                # neutral beam power deposition
            'power_NB':['S', 'W', ['power_nbi'] ],
            'ni_NB(rho)':['P', 'm^-3', ['nbeami'],['~rho_nbi'] ],
            'Pe_NB_dens(rho)':['P', 'W m^-3', ['pbe','rho_nbi', 'rho_eq', 'vol'],['rho_nbi'] ],
            'Pi_NB_dens(rho)':['P', 'W m^-3', ['pbi','rho_nbi', 'rho_eq', 'vol'],['rho_nbi'] ],
            'Pe_NB(rho)':['P', 'W', ['pbe'],['~rho_nbi'] ],
            'Pi_NB(rho)':['P', 'W', ['pbi'],['~rho_nbi'] ],
            'Pe_NB_cum(rho)':['P', 'W', ['pbe'],['~rho_nbi'] ],
            'Pi_NB_cum(rho)':['P', 'W', ['pbi'],['~rho_nbi'] ],
            'Pe_NB_total':['S', 'W', ['pbe'] ],
            'Pi_NB_total':['S', 'W', ['pbi'] ],
            'J_NB(rho)':['P', 'A m^-2', ['curbeam','rho_nbi', 'rho_eq', 'area'],['rho_nbi'] ],
            'I_NB(rho)':['P', 'A', ['curbeam'],['~rho_nbi'] ],
            'I_NB_cum(rho)':['P', 'A', ['curbeam'],['~rho_nbi'] ],
            'I_NB_total':['S', 'A', ['curbeam'] ],
            'Eperp_NBI(rho)':['P', 'keV', ['eperp_beami'],['~rho_nbi'] ],
            'Epll_NBI(rho)':['P', 'keV', ['epll_beami'],['~rho_nbi'] ]
}

    
    # Dictionary of dependencies that whose names appear in plasma state lists rather than
    # at the top level of the plasma state. The prime example is a thermal ion species whose
    # names appear in the plasma state list S_name.
    
labeled_ps_quantities = {
            'n_H(0)':['S_name', 'H'],
            'n_H_ave':['S_name', 'H'],
            'n_H(rho)':['S_name', 'H'],
            'n_D(0)':['S_name', 'D'],
            'n_D_ave':['S_name', 'D'],
            'n_D(rho)':['S_name', 'D'],
            'n_T(0)':['S_name', 'T'],
            'n_T_ave':['S_name', 'T'],
            'n_T(rho)':['S_name', 'T'],
            'n_He4(0)':['S_name', 'He4'],
            'n_He4_ave':['S_name', 'He4'],
            'n_He4(rho)':['S_name', 'He4']
            }

if debug:
    monitorDefinition.update( {    # testing
            'dummy':['S', 'arb', ['stuff']], 'q':['P', 'arb',['ns'],['rho'] ]
    } )

print 'monitor_comp_version = ', monitor_comp_version
print 'metaData = ',monitorDefinition['monitor_comp_metaData']


class monitor(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)


#----------------------------------------------------------------------------------------------
#
# Define some utility functions
#
#----------------------------------------------------------------------------------------------
    
    def find_label_index_in_list(self, varDepsDict, list_name, label):  
    
        '''
        list_name = name of a plasma state varaible that is an item list e.g. species name list
        label = string = a name that may appear in the list
        Returns the index of label in the list if it's there.  Raises exception otherwise.
        '''
        ps_list = varDepsDict[list_name][:]
        str_list = map(lambda x: ''.join(x).strip(), ps_list)
        try:
            i = str_list.index(label)
        except Exception, e:
            print 'did not find item', label, ' in list ', list_name
            raise Exception
        return i

#----------------------------------------------------------------------------------------------
#
# Define some functions used in calculating monitor variables
#
#----------------------------------------------------------------------------------------------

    def linear_grid_interp (self, y, g1, g2):
    
        '''
        Simple Q&D linear interpolator for profiles, y, defined on Plasma State rho grid,
        g1, to rho grid, g2.  Both grids must begin at 0.0 and end at 1.0 and be 
        monotonically increasing.
        ''' 
        # Check validity of inputs
        if g2[0] != 0.0 or g1[0] != 0.0 or g2[-1] != 1.0 or  g1[-1] != 1.0 :
            print 'Improper range for grids'
            print 'grid1[0] = ', g1[0], '  grid1[-1] = ', g1[-1]
            print 'grid2[0] = ', g2[0], '  grid2[-1] = ', g2[-1]
            return 1
        if len(y) != len(g1):
            print 'linear_grid_interp Error: length of y = ', len(y), \
            ' does not match length of g1 = ', len(g1)
            return 1
    
        yi = []
        yi.append(y[0])
    
        j_up = 1
        for i in range(1, len(g2)-1):
            while g2[i] > g1[j_up]:  # Find the index of the point on the g1 grid just above the 
                j_up = j_up + 1      # g2 point
    
            interp = y[j_up-1] + (g2[i] - g1[j_up-1])*(y[j_up] - y[j_up-1])/\
            (g1[j_up] - g1[j_up-1])
            yi.append(interp)
        yi.append(y[-1])
    
        return yi
#--------------------------------------------------------------------------------------------

    def area_grid_interp (self, y, g1, g2):
    
        '''
        Simple Q&D area interpolator for profiles, y, defined on Plasma State rho grid, g1,
        g1, to rho grid, g2.  Both grids must begin at 0.0 and end at 1.0 and be monotonically
        increasing.  This also works for volume since in a torus volume goes like 2*Pi*area.  
        Just use vol for argument y.
        
        There is a problem when linearly interpolating a quantity which increases almost as 
        rho^2 onto a finer grid (like ps%area = area enclosed in a flux surface). To calculate 
        a density (m^-2) from a zone based variable, like current/zone, one needs to divide
        by area/zone.  In a linear interpolation onto the finer grid the area/zone is the 
        same in all the sub-zones and takes a sudden jump when crossing over a grid point of 
        the original coarse grid.  In fact the area/zone should increase linearly in the  
        sub-zones to reflect the quadratic profile.  This routine fits a function of the from
        A + B*rho^2 to each pair of points on the original grid and interpolates that.  It 
        might be better to fit a quadratic to each 3 points, but that complicates the ends 
        and doesn't seem worth it now.
        ''' 
        # Check validity of inputs
                
        if abs(g1[0]) < 1.e-12:
        	g1[0] = 0.
        if abs(g2[0]) < 1.e-12:
        	g2[0] = 0.
        if abs(g1[-1] -1.0) < 1.e-12:
        	g1[-1] = 1.0
        if abs(g2[-1] -1.0) < 1.e-12:
        	g2[-1] = 1.0
        
        if g2[0] != 0.0 or g1[0] != 0.0 or g2[-1] != 1.0 or  g1[-1] != 1.0 :
            print 'Improper range for grids'
            print 'grid1[0] = ', g1[0], '  grid1[-1] = ', g1[-1]
            print 'grid2[0] = ', g2[0], '  grid2[-1] = ', g2[-1]
            return 1
        if len(y) != len(g1):
            print 'linear_grid_interp Error: length of y = ', len(y), \
            ' does not match length of g1 = ', len(g1)
            return 1
    
        yi = []
        yi.append(y[0])
    
        j_up = 1
        for i in range(1, len(g2)-1):
            while g2[i] > g1[j_up]:  # Find the index of the point on the g1 grid just above the 
                j_up = j_up + 1      # g2 point
    
            interp = y[j_up-1] + (g2[i]**2 - g1[j_up-1]**2)*(y[j_up] - y[j_up-1])/\
            (g1[j_up]**2 - g1[j_up-1]**2)
            yi.append(interp)
        yi.append(y[-1])
    
        return yi

#--------------------------------------------------------------------------------------------


    def vol_average(self, x, rho_x, vol, rho_eq):

        # Interpolate vol(rho_eq) onto rho grid
        vol_rho_x = self.area_grid_interp(vol, rho_eq, rho_x)

        vol_ave = 0.0
        for i in range(len(x)-1):
            vol_ave = vol_ave + (vol_rho_x[i+1] - vol_rho_x[i])*x[i]
        vol_ave = vol_ave/vol_rho_x[-1]
        return vol_ave

#--------------------------------------------------------------------------------------------

    def cumulate(self, x):
    
        c = x[0]
        cum = [c]
        for i in range(1,len(x)):
            c = c + x[i]
            cum.append(c)
        return cum

#--------------------------------------------------------------------------------------------

    def zone_to_surf_dens(self, x, rho_x, area, rho_eq):
        # Takes a zone based quantity, x, defined on grid, rho_x, and converts it
        # to a density/area.  Density at rho_x = 0 is taken as density in the first
        # zone.  Density at rho_x = 1 is taken as density in last zone.  At interior
        # grid points there is a jump in density, so take the mean of the zone above
        # and the zone below.  This also works for volume densities since in a torus
        # volume goes like 2*Pi*area.  Just use vol for area argument.
        
        # NB: Remember the dimension of zone based quantities is 1 less than the
        #     number of grid points
        
        # Interpolate area(rho_eq) onto rho_x grid
        area_rho_x = self.area_grid_interp(area, rho_eq, rho_x)
        
        dens = [ x[0]/area_rho_x[1] ]
        for i in range(len(x)-1):
            dens_m = x[i]/(area_rho_x[i+1] - area_rho_x[i])
            dens_p = x[i+1]/(area_rho_x[i+2] - area_rho_x[i+1])
            dens.append( 0.5*(dens_m + dens_p) )
        dens.append( x[-1]/(area_rho_x[-1] - area_rho_x[-2]) )
        return dens

#--------------------------------------------------------------------------------------------

    def line_average_density(self, dens, rho, R_surfMax, R_surfMin, rho_eq_geo):
    
        # Interpolate R_surfMax(rho_eq_geo) and R_surfMin(rho_eq_geo) onto rho grid
        R_Max_rho = self.linear_grid_interp(R_surfMax, rho_eq_geo, rho)
        R_Min_rho = self.linear_grid_interp(R_surfMin, rho_eq_geo, rho)
    
        line_ave = 0.0
        for i in range(len(dens)-1):
            line_ave = line_ave + (R_Max_rho[i+1] - R_Max_rho[i] - R_Min_rho[i+1]\
                       + R_Min_rho[i])*dens[i]
        line_ave = line_ave/(R_Max_rho[-1] - R_Min_rho[-1])
        
        if debug:
            print 'line_ave = ', line_ave
        return line_ave

#--------------------------------------------------------------------------------------------

    def calculate_qPsi(self, Psi, q_eq, rho_eq, psipol):

        #Find rho_eq for psi_poloidal(rho)/ psi_poloidal(1) = Psi, where 0 < Psi < 1
        i_low = len(psipol)- 1
        while psipol[i_low]/psipol[-1] > Psi:
            i_low = i_low - 1
        rho_l = rho_eq[i_low]
        rho_u = rho_eq[i_low + 1]
        psi_l = psipol[i_low]
        psi_u = psipol[i_low +1]
        rho_eq_Psi = rho_l + (rho_u - rho_l)*( 0.95*psipol[-1] - psi_l)/(psi_u - psi_l)

        #interpolte to get q_Psi
        q_Psi = q_eq[i_low] + (rho_eq_Psi - rho_l)*(q_eq[i_low+1] -q_eq[i_low])/(rho_u - rho_l)
        return q_Psi

#--------------------------------------------------------------------------------------------

    def calculate_triangPsi(self, Psi, triang, psipol, rho_eq, rho_eq_geo):

        #Find rho_eq for psi_poloidal(rho)/ psi_poloidal(1) = Psi, where 0 < Psi < 1
        i_low = len(psipol)- 1
        while psipol[i_low]/psipol[-1] > Psi:
            i_low = i_low - 1
        rho_l = rho_eq[i_low]
        rho_u = rho_eq[i_low + 1]
        psi_l = psipol[i_low]
        psi_u = psipol[i_low +1]
        rho_eq_Psi = rho_l + (rho_u - rho_l)*( Psi*psipol[-1] - psi_l)/(psi_u - psi_l)

        #Find rho_eq_geo for psi_poloidal(rho)/ psi_poloidal(1) = Psi, where 0 < Psi < 1
        i_low = len(rho_eq_geo)- 1
        while rho_eq_geo[i_low] > rho_eq_Psi:
            i_low = i_low - 1
        rho_l = rho_eq_geo[i_low]
        rho_u = rho_eq_geo[i_low + 1]
        triang_l = triang[i_low]
        triang_u = triang[i_low +1]

        #interpolte to get triang_Psi
        triang_Psi = triang_l + (rho_eq_Psi - rho_l)*(triang_u -triang_l)/(rho_u - rho_l)
        return triang_Psi

#--------------------------------------------------------------------------------------------

    def calculate_elongPsi(self, Psi, elong, psipol, rho_eq, rho_eq_geo):

        # Find rho_eq for psi_poloidal(rho)/ psi_poloidal(1) = Psi, where 0 < Psi < 1
        i_low = len(psipol)- 1
        while psipol[i_low]/psipol[-1] > Psi:
            i_low = i_low - 1
        rho_l = rho_eq[i_low]
        rho_u = rho_eq[i_low + 1]
        psi_l = psipol[i_low]
        psi_u = psipol[i_low +1]
        rho_eq_Psi = rho_l + (rho_u - rho_l)*( Psi*psipol[-1] - psi_l)/(psi_u - psi_l)

        # Find rho_eq_geo for psi_poloidal(rho)/ psi_poloidal(1) = Psi, where 0 < Psi < 1
        i_low = len(rho_eq_geo)- 1
        while rho_eq_geo[i_low] > rho_eq_Psi:
            i_low = i_low - 1
        rho_l = rho_eq_geo[i_low]
        rho_u = rho_eq_geo[i_low + 1]
        elong_l = elong[i_low]
        elong_u = elong[i_low +1]

        # interpolte to get elong_Psi
        elong_Psi = elong_l + (rho_eq_Psi - rho_l)*(elong_u -elong_l)/(rho_u - rho_l)
        return elong_Psi

#--------------------------------------------------------------------------------------------

    def calculate_li_3(self, curt, q_eq, grho2r2i, phit, R_axis, vol):

        mu0 = 4.0*math.pi/10000000.0

        # Do volume integral - N.B. Everything in the integral is defined on the
        # rho_eq grid except rho.  So use rho**2 = Phi_toroidal/Phi_toroidal_edge, then
        # everything is on rho_eq grid.

        li_3 = 0.0
        for i in range(len(q_eq)-1):
            li_3 = li_3 + ( vol[i+1] - vol[i] )*phit[i]*grho2r2i[i]/(q_eq[i]*q_eq[i])
        li_3 = 8.0*li_3/(R_axis * mu0*mu0 * curt[-1]*curt[-1])

        if debug:
            print 'li_3 = ', li_3
        return li_3

#--------------------------------------------------------------------------------------------

    def calculate_N_GW(self, ns,rho, rho_eq_geo, R_surfMin, R_surfMax, curt):
    
        ne_rho = ns[0]
        line_ave = self.line_average_density(ne_rho, rho, R_surfMax, R_surfMin, rho_eq_geo)
        # normalize ne to units of 10**20 and Ip to MA
        line_ave = line_ave/1.0e20
        Ip = curt[-1]/1.0e6
        a = (R_surfMax[-1] - R_surfMin[-1])/2.0
        N_GW = math.pi*a**2*line_ave/Ip
        
        if debug:
            print 'N_GW = ', N_GW
        return N_GW
    


# ------------------------------------------------------------------------------
#
# Stages the input files (currently there are none) for the monitor executable  and
# inital plasma state file into monitor component work directory.  Calls the
# monitor executable in "INIT" mode, which generates the monitor file
# with the initial data.  Copies 'monitor_file.nc' generated by the executable into
# the SWIM w3 directory (/p/swim/w3_html) with the runid prepended to the file name
#
# ------------------------------------------------------------------------------

    def init(self, timeStamp=0):
        """Function init generates the initial netcdf file for monitor data
        """

        print 'monitor_comp.init() called'

        services = self.services

        time.sleep(5)
        self.run_id = services.get_config_param('PORTAL_RUNID')
        #self.run_id = self.get_config_param(services,'PORTAL_RUNID')
        print 'run_id = ', self.run_id
    	print 'monitor file = ', monitor_fileName
    	print 'monitor pdf file = ', pdf_fileName

        self.cdfFile = self.run_id+'_' + monitor_fileName
        self.pdfFile = self.run_id+'_' + pdf_fileName
        services.log('w3 monitor file = ' + self.cdfFile)
        services.log('state pdf file = ' + self.pdfFile)

        BIN_PATH = self.get_component_param(services, 'BIN_PATH')
        self.PCMF_bin = os.path.join(self.BIN_PATH, 'PCMF.py')
         
    # Copy current state over to working directory
        services.stage_plasma_state()

        services.stage_input_files(self.INPUT_FILES)
        cur_state_file = services.get_config_param('CURRENT_STATE')

    # Generate initial monitor file
        retcode = self.init_monitor_file(cur_state_file, timeStamp)
        if (retcode != 0):
            services.log('Error executing INIT:  init_monitor_file')
            return 1

    # copy monitor file to w3 directory
        try:
            shutil.copyfile(monitor_fileName,
                            os.path.join(self.W3_DIR, self.cdfFile))
        except IOError, (errno, strerror):
            print 'Error copying file %s to %s: %s' % \
                (monitor_file, self.cdfFile, strerror)

   # Generate pdf file with PCMF.py
        cmd = [self.PCMF_bin, monitor_fileName]
        print 'Executing = ', cmd
        services.send_portal_event(event_type = 'COMPONENT_EVENT',\
          event_comment =  cmd)
        retcode = subprocess.call(cmd)
        if (retcode != 0):
            logMsg = 'Error executing '.join(map(str, cmd))
            self.services.error(logMsg)
            raise Exception(logMsg)

   # copy pdf file to w3 directory
        try:
            shutil.copyfile(pdf_fileName,
                            os.path.join(self.W3_DIR, self.pdfFile))
        except IOError, (errno, strerror):
            print 'Error copying file %s to %s: %s' % \
                (pdf_fileName, self.pdfFile, strerror)

    # Copy config file to w3 directory
        conf_file = services.get_config_param('SIMULATION_CONFIG_FILE')
        print 'conf_file = ', conf_file
        conf_file_name = os.path.split(conf_file)[1]
        new_file_name = self.run_id + '_' + conf_file_name
        new_full_path = os.path.join(self.W3_DIR, new_file_name)
        try:
            shutil.copyfile(conf_file, new_full_path)
        except IOError, (errno, strerror):
            print 'Error copying file %s to %s: %s' % \
                (conf_file, new_full_path, strerror)

        return

# ------------------------------------------------------------------------------
#
# restart function
#
# Load the internal state needed to restart the monitor component.  In particular
# monitorVars and ps_VarsList are unpickeld pickled from a file "monitor_restart".
#
# ------------------------------------------------------------------------------

    def restart(self, timeStamp):
        """
        Function restart loads the internal monitor state data needed for
        restart_MonitorComponent
        """
        print 'monitor_comp.restart() called'

        services = self.services
        global monitorVars, ps_VarsList, monitorDefinition
        
        self.run_id = self.get_config_param(services,'PORTAL_RUNID')

        self.cdfFile = self.run_id+'_' + monitor_fileName
        services.log('w3 monitor file = ' + self.cdfFile)
        
    # Get restart files listed in config file.        
        try:
            restart_root = services.get_config_param('RESTART_ROOT')
            restart_time = services.get_config_param('RESTART_TIME')
            services.get_restart_files(restart_root, restart_time, self.RESTART_FILES)
        except Exception, e:
            print 'Error in call to get_restart_files()' , e
            raise

    # copy monitor file to w3 directory
        try:
            shutil.copyfile(monitor_fileName,
                            os.path.join(self.W3_DIR, self.cdfFile))
        except IOError, (errno, strerror):
            print 'Error copying file %s to %s: %s' % \
                (monitor_fileName, self.cdfFile, strerror)

   # Generate pdf file with PCMF.py
        cmd = [self.PCMF_bin, monitor_fileName]
        print 'Executing = ', cmd
        services.send_portal_event(event_type = 'COMPONENT_EVENT',\
          event_comment =  cmd)
        retcode = subprocess.call(cmd)
        if (retcode != 0):
            logMsg = 'Error executing '.join(map(str, cmd))
            self.services.error(logMsg)
            raise Exception(logMsg)

    # copy pdf file to w3 directory
        try:
            shutil.copyfile(pdf_fileName,
                            os.path.join(self.W3_DIR, self.pdfFile))
        except IOError, (errno, strerror):
            print 'Error copying file %s to %s: %s' % \
                (pdf_fileName, self.pdfFile, strerror)
    
    # Load monitorVars, monitorDefinition and ps_VarsList from pickle file "monitor_restart".

        pickleDict = {'monitorVars' : monitorVars, 'ps_VarsList': ps_VarsList,\
                     'monitorDefinition':monitorDefinition}
        pickFile = open('monitor_restart', 'r')
        pickleDict = pickle.load(pickFile)
        pickFile.close()
        monitorVars = pickleDict['monitorVars']
        ps_VarsList = pickleDict['ps_VarsList']
        monitorDefinition = pickleDict['monitorDefinition']
        print 'monitorDefinition = ', monitorDefinition
        
        print 'monitor restart finished'
        return 0

# ------------------------------------------------------------------------------
#
# step function
#
# Stages current input files (primarily plasma state). Updates the monitor_file.nc
# from current plasma state.  And copies updated monitor file to w3 directory
#
# ------------------------------------------------------------------------------

    def step(self, timeStamp):
        print '\nmonitor_comp.step() called'

        services = self.services

        if (self.services == None) :
            print 'Error in monitor_comp: step() : no framework services'
            return 1

#        monitor_file = 'monitor_file.nc'

    # Copy current and prior state over to working directory
        services.stage_plasma_state()
        cur_state_file = services.get_config_param('CURRENT_STATE')

        monitor = os.path.join(self.BIN_PATH, 'monitor_comp')

    # Call Load new data into monitor file
        retcode = self.update_monitor_file(cur_state_file, timeStamp)
        if (retcode != 0):
            print 'Error executing command: ', monitor
            return 1

    # "Archive" output files in history directory
        services.stage_output_files(timeStamp, self.OUTPUT_FILES)

    # copy monitor file to w3 directory
        try:
            shutil.copyfile(monitor_fileName,
                            os.path.join(self.W3_DIR, self.cdfFile))
        except IOError, (errno, strerror):
            print 'Error copying file %s to %s: %s' % \
                (monitor_fileName, self.cdfFile, strerror)

   # Generate pdf file with PCMF.py
        cmd = [self.PCMF_bin, monitor_fileName]
        print 'Executing = ', cmd
        services.send_portal_event(event_type = 'COMPONENT_EVENT',\
          event_comment =  cmd)
        retcode = subprocess.call(cmd)
        if (retcode != 0):
            logMsg = 'Error executing '.join(map(str, cmd))
            self.services.error(logMsg)
            raise Exception(logMsg)

    # copy pdf file to w3 directory
        try:
            shutil.copyfile(pdf_fileName,
                            os.path.join(self.W3_DIR, self.pdfFile))
        except IOError, (errno, strerror):
            print 'Error copying file %s to %s: %s' % \
                (pdf_fileName, self.pdfFile, strerror)

# ------------------------------------------------------------------------------
#
# checkpoint function
#
# Saves restart files to restart directory.  Should include: monitor_restart and
# monitor_file.nc
# ------------------------------------------------------------------------------

    def checkpoint(self, timestamp=0.0):
        print 'monitor.checkpoint() called'
        services = self.services
        services.save_restart_files(timestamp, self.RESTART_FILES)
        

# ------------------------------------------------------------------------------
#
# finalize function
#
# Calls monitor executable in "FINALIZE" mode which sets
# 'running' atribute to false so Elvis will stop watching the monitor file.
#
# ------------------------------------------------------------------------------

    def finalize(self, timestamp=0.0):

        # Open the monitor file for output
        monitor_file = Dataset(monitor_fileName, 'r+', format = 'NETCDF3_CLASSIC')

        # set 'running' attribute to tell ELVis to stop watching this file
        setattr(monitor_file, 'running', 'false')

        # Close plasma_state
        monitor_file.close()

        print 'monitor finalize finished'
        return 0


#----------------------------------------------------------------------------------------------
#
# init_monitor_file function
#
# Analyze the Plasma State dependencies of the variables requested to be monitored.
# Define the initial monitor netcdf file.  Load the initial data into the monitor file
# By a call to the 'step' function.  Also saves the internal state needed to restart
# the monitor component.  In particular monitorVars and ps_VarsList are pickled to a
# file called monitor_restart.
#----------------------------------------------------------------------------------------------

    def init_monitor_file(self, cur_state_file, timeStamp = 0):
        """Function init_monitor_file generates the initial netcdf file for monitor data
           and saves the internal state needed to restart the monitor component
        """
        global requestedVars, monitorVars, monitorDims, ps_VarsList, monitorGrids,\
        density_based_profiles, zone_based_profiles, cumulative_profiles,\
        include_density_based_profiles, include_zone_based_profiles,\
        include_cumulative_profiles,\
        monitorDefinition

        print ' '
        print 'monitor_component: init_monitor_file'
        
        # Check if include_density_based_profiles, include_zone_based_profiles, and
        # include_cumulative_profiles are set with config variables.  Otherwise use
        # defaults defined above.

        try:
            include_density_based_profiles = self.INCLUDE_DENSITY_BASED_PROFILES
        except:
            pass
        try:
            include_zone_based_profiles = self.INCLUDE_ZONE_BASED_PROFILES
        except:
            pass
        try:
            include_cumulative_profiles = self.INCLUDE_CLUMULATIVE_PROFILES
        except:
            pass

        if include_density_based_profiles:
            print 'include_density_based_profiles'
            requestedVars = requestedVars + density_based_profiles
        if include_zone_based_profiles:
            print 'include_zone_based_profiles'
            requestedVars = requestedVars + zone_based_profiles
        if include_cumulative_profiles:
            print 'include_cumulative_profiles'
            requestedVars = requestedVars + cumulative_profiles

        if debug:
            requestedVars = requestedVars + ['dummy', 'crazy']

        print 'cur_state_file =', cur_state_file
        print 'timeStamp =', timeStamp
        print 'requestedVars = ', requestedVars

        #Open Plasma State file
        plasma_state = Dataset(cur_state_file, 'r', format = 'NETCDF3_CLASSIC')
    
        all_ps_DimNames = plasma_state.dimensions.keys()
    
        all_ps_Dims = plasma_state.dimensions
    
        all_ps_VarNames = plasma_state.variables.keys()
    
        Global_label = str(chartostring(plasma_state.variables['Global_label'][:]))
        print 'Global_label = ', Global_label
    
        RunID = str(chartostring(plasma_state.variables['RunID'][:]))
        print 'RunID = ', RunID
    
        tokamak_id = str(chartostring(plasma_state.variables['tokamak_id'][:]))
        print 'tokamak_id = ', tokamak_id
    
        shot_number = plasma_state.variables['shot_number'].getValue()

        if debug:
            print ' '
            all_ps_DimNames.sort()
            print 'all_ps_Dims = ', all_ps_DimNames
            print ' '
            all_ps_VarNames.sort()
            print 'all_ps_VarNames = ', all_ps_VarNames
            print ' '

        # Some monitor variables, e.g Q_nuc_FUS, may make sense with or without some of the
        # plasma state variables present.  For example Q_nuc_FUS is well defined even if some
        # of the heating sources (e.g. ECH) are not present.  So modify the dependency list
        # in the monitorDefinition dictionary for these variables based on what is available
        # in the plasma state.
        
        # Q_NUC_FUS
        if 'power_ec' in all_ps_VarNames:
            monitorDefinition['Q_nuc_FUS'][2].append('power_ec')
        if 'power_lh' in all_ps_VarNames:
            monitorDefinition['Q_nuc_FUS'][2].append('power_lh')
        if 'power_ic' in all_ps_VarNames:
            monitorDefinition['Q_nuc_FUS'][2].append('power_ic')
        if 'power_nbi' in all_ps_VarNames:
            monitorDefinition['Q_nuc_FUS'][2].append('power_nbi')

        #
        ## Determine if data is available to calculate all requested monitor variables
        #
        for var in requestedVars:
            var_deps_ok = True

            # Check if var is in monitorDefinition .i.e we know how to calculate it
            if var not in monitorDefinition.keys():
                print ' '
                print 'Requested variable ', var, \
                      ' not in monitorDefinition.keyes.'
                var_deps_ok = False
                continue

            if debug:
                print ' '
                print 'requestedVar = ', var
                print ' '
                print var, 'monitorDefinition = ', monitorDefinition[var]

            varKind = monitorDefinition[var][0]
            varDepsList = monitorDefinition[var][2]

            newVars = []
            newGrids = []
            newDims = []
            for dep in varDepsList:
                # Check if needed variable is in Plasma State
                if dep in all_ps_VarNames:
                    # Check if already in full list
                    if dep not in ps_VarsList:
                        newVars.append(dep)
                    # Even if the dependency is in the Plasma State it might a labeled list.  
                    # Check if var depends on a labeled list and if this dep is that labeled list.
                    if var in labeled_ps_quantities.keys():
                        if dep == labeled_ps_quantities[var][0]:
                            # If this is a labeled list that var depends on check if all the items
                            # that var depends on are in the labeled list.
                            ps_list = plasma_state.variables[dep][:]
                            str_list = map(lambda x: ''.join(x).strip(), ps_list)
                            if debug:
                                print 'list variable = ', var, '  str_list = ', str_list
                            for label in labeled_ps_quantities[var][1:]:
                                if label not in str_list :
                                    print 'did not find label ', label, ' in Plasma State list ', dep
                                    var_deps_ok = False
                else:
                    print ' '
                    print var, 'dependency ', dep, 'not found in Plasma State'
                    var_deps_ok = False
                    break

            # for profiles check that needed grids are in plasma state
            if varKind in ['P','2D']:

                if debug:
                    print var, ' kind = ', varKind

                varGridsList = monitorDefinition[var][3]

                # For profiles and 2D check for right number of grids present

                if varKind == 'P' and len(varGridsList) != 1:
                    var_deps_ok = False
                    print ' '
                    print len(varGridsList), ' grids specified for ', var, ' should = 1'
                    break

                if varKind == '2D' and len(varGridsList) != 2:  # check for right number of grids
                    var_deps_ok = False
                    print ' '
                    print len(varGridsList), ' grids specified for ', var, ' should = 2'
                    break

                if debug:
                    print ' '
                    print 'varGridsList = ', varGridsList

                # Check if grid dependencies are present
                for grid in varGridsList:
                    # Delete the tilde, if any, so can find in plasma state
                    if grid[0] == '~':
                        ps_grid_name = grid[1:]
                    else:
                        ps_grid_name = grid
                    # Check if needed grid is in Plasma State
                    if ps_grid_name in all_ps_VarNames:
                        # Check if already in full list, if not append it to list
                        if grid not in monitorGrids:
                            newGrids.append(grid)
                    else:
                        print ' '
                        print var, 'grid ', ps_grid_name, \
                              'not found in Plasma State'
                        var_deps_ok = False
                        break

            # for Scalar Lists check that needed dimensions are in plasma state
            if varKind in ['SL']:

                if debug:
                    print var, ' kind = ', varKind

                varDimList = monitorDefinition[var][3]

                # For scalar lists check for right number of dimensions = 1

                if len(varDimList) != 1:
                    var_deps_ok = False
                    print ' '
                    print len(varDimList), ' dimensons specified for ', var, ' should = 1'
                    break

                if debug:
                    print ' '
                    print 'varDimList = ', varDimList

                # Check if dim dependencies are present
                for dim_name in varDimList:
                    # Check if needed grid is in Plasma State
                    if dim_name in all_ps_DimNames:
                        # Check if already in full list, if not append it to list
                        if dim_name not in monitorDims:
                            newDims.append(dim_name)
                    else:
                        print ' '
                        print var, 'dim ', dim_name, \
                              'not found in Plasma State'
                        var_deps_ok = False
                        break

            # Keep this monitor variable if all dependencies are ok
            if var_deps_ok == True:
                monitorVars = monitorVars + [var]
                # add to list of needed ps variables and grids
                ps_VarsList = ps_VarsList + newVars
                monitorGrids = monitorGrids + newGrids
                monitorDims = monitorDims + newDims

        print ' '
        print '\nmonVars = ', monitorVars
        print '\nps_VarsList = ', ps_VarsList
        print '\nmonitorGrids = ', monitorGrids
        print '\nmonitorDims = ', monitorDims



        #
        ## Define monitor file
        #

        # Open the monitor file for output
        monitor_file = Dataset(monitor_fileName, 'w', format = 'NETCDF3_CLASSIC')

        # make global netcdf attribute of monitor component version number
        setattr(monitor_file, 'monitor_comp_version', monitor_comp_version)

        # make global netcdf attribute to tell ELVis to watch this file
        setattr(monitor_file, 'running', 'true')

        # make global netcdf attribute for Global_label
        setattr(monitor_file, 'Global_label', Global_label)

        # make global netcdf attribute for RunID
        setattr(monitor_file, 'RunID', RunID)

        # make global netcdf attribute for tokamak_id
        setattr(monitor_file, 'tokamak_id', tokamak_id)

        # make global netcdf attribute for shot_number
        setattr(monitor_file, 'shot_number', shot_number)

        # Create unlimited time dimension and define time variable
        monitor_file.createDimension('timeDim', None)
        time = monitor_file.createVariable('time', float, ('timeDim',))
        setattr(time, 'units', 'sec')

        # Create grid dimensions and variables and load up grids
        # use this to keep track of mon_grid name, grid_dim_name and grid_dim_value
        grid_map = {}
        for grid in monitorGrids:
            # determine if this is a zone based grid (i.e. length in ps -1)
            # Delete the tilde, if any, so can find in plasma state
            if grid[0] == '~':
                ps_grid_name = grid[1:]
                mon_grid_name = ps_grid_name + '_zone'
            else:
                ps_grid_name = grid
                mon_grid_name = ps_grid_name

            # Get the netcdf variable object for this grid
            grid_obj = plasma_state.variables[ps_grid_name]
            grid_dim_name = grid_obj.dimensions
            grid_dim_value = len(plasma_state.dimensions[grid_dim_name[0]])
            grid_val = grid_obj[:]
            grid_units = getattr(grid_obj, 'units')
            grid_shape = grid_obj.shape

            if debug:
                print ' '
                print ps_grid_name, 'dim_name = ', grid_dim_name
                print ps_grid_name, 'dim_value = ', grid_dim_value
                print ps_grid_name, 'value = ', grid_val
                print ps_grid_name, 'units = ', grid_units
                print ps_grid_name, 'shape = ', grid_shape

            if len(grid_shape) != 1:
                print 'whoa!  grid ', ps_grid_name, \
                      ' shape not 1D.  I don''t know what to do'
                return 1

            # case of zone boundary grid (copy straight from plasma state)
            if mon_grid_name == ps_grid_name:
                monitor_file.createDimension(grid_dim_name[0],
                                             grid_dim_value)
                mon_obj = monitor_file.createVariable(mon_grid_name,
                                                      float,
                                                      grid_dim_name )
                mon_obj[:] = grid_val

            else: # case of zone centered grid (average grid boundaries)
                grid_dim_name = ('dm1'+grid_dim_name[0][3:],)
                grid_dim_value = grid_dim_value - 1
                # print 'modified grid_dim_name = ', grid_dim_name
                monitor_file.createDimension(grid_dim_name[0], grid_dim_value)
                mon_obj = monitor_file.createVariable(mon_grid_name,
                                                      float,
                                                      grid_dim_name )
                mon_obj[:] = 0.5 * ( grid_val[0: grid_dim_value] + \
                                     grid_val[1: grid_dim_value + 1] )

            #Generate units atribute for grid
            setattr(mon_obj, 'units', grid_units)
            grid_map[grid] = [mon_grid_name, grid_dim_name, grid_dim_value]

        # Create other non-grid dimensions e.g. for scalar lists
        for dim_name in monitorDims:

            dim_value = len(plasma_state.dimensions[dim_name])
            monitor_file.createDimension(dim_name, dim_value)

        # Create monitor variables
        for var in monitorVars:

            if debug:
                print 'creating variable ', var

            # Generate the dimension tuple
            dims = ('timeDim',)
            varKind = monitorDefinition[var][0]

            if varKind == 'S':
                mon_obj = monitor_file.createVariable(var, float, dims)

            elif varKind == 'P':
                varGridsList = monitorDefinition[var][3]
                # add in dimension names for plotting grids
                for grid in varGridsList:
                    dims = dims + grid_map[grid][1]

                mon_obj = monitor_file.createVariable(var, float, dims )

            elif varKind == 'SL':
                varDimList = monitorDefinition[var][3]
                # add in dimension names
                for dim_name in varDimList:
                    dims = dims + (dim_name,)

                mon_obj = monitor_file.createVariable(var, float, dims )

            elif varKind == '2D':   # See note at top about the order of dims in definition
                varGridsList = monitorDefinition[var][3]
                for grid in varGridsList:  # add in dimension names for plotting grids
                    dims = dims + grid_map[grid][1]

                if debug:
                    print var, ' dimensions = ', dims
                mon_obj = monitor_file.createVariable(var, float, dims )

            varUnits = monitorDefinition[var][1]
            #Generate units atribute for var
            setattr(mon_obj, 'units', varUnits)

        # Finished defining monitor file.  Close it.
        monitor_file.close()

        # Close plasma_state
        plasma_state.close()

        # insert intitial data
        self.step(timeStamp)

        print 'monitor file initialization finished'
        
    # Save monitorVars and ps_VarsList and monitorDefinition, are pickled to file "monitor_restart".
    
        pickleDict = {'monitorVars' : monitorVars, 'ps_VarsList': ps_VarsList,\
                     'monitorDefinition':monitorDefinition}
        pickFile = open('monitor_restart', 'w')
        pickle.dump(pickleDict, pickFile)
        pickFile.close() 

        return 0

#----------------------------------------------------------------------------------------------
#
# update_monitor_file function
#
# Opens current plasma state file and monitor netcdf file.  Pulls the needed variables out of
# the plasma state file.  Computes the monitor variables from plasma state data (in functinon
# calculate_MonitorVariable).  And loads the data into the monitor netcdf variables.
#
#----------------------------------------------------------------------------------------------

    def update_monitor_file(self, cur_state_file, timeStamp = 0):

        print ' '
        print 'monitor_component: update_monitor_file'

        #Open Plasma State file
        plasma_state = Dataset(cur_state_file, 'r', format = 'NETCDF3_CLASSIC')

        # Get all Plasma State variable objects for dependencies
        ps_variables = {}
        for var in ps_VarsList:
            ps_variables[var] = plasma_state.variables[var]

        if debug:
            print ' '
            print 'ps_variables.keys() = ', ps_variables.keys()

        # Open the monitor file for output
        monitor_file = Dataset(monitor_fileName, 'r+', format = 'NETCDF3_CLASSIC')

        if debug:
            all_mon_Dims = monitor_file.dimensions
            all_mon_VarNames = monitor_file.variables.keys()
            print ' '
            print 'all_mon_Dims = ', all_mon_Dims
            print ' '
            print 'all_mon_VarNames = ', all_mon_VarNames

        # Get time variable object
        time = monitor_file.variables['time']
        n_step =time.shape[0]    # Time step number (initially 0)
        print 'time step number = ', n_step
        time[n_step] = float(timeStamp)

        # Insert data into monitor variables
        for var in monitorVars:
            # Get the netcdf variable object for this variable
            var_obj = monitor_file.variables[var]

            # Get Plasma State variables for this monitor variable's dependencies
            varDepsDict ={}
            varDepsList = monitorDefinition[var][2]
            for dep in varDepsList:
                varDepsDict[dep] = ps_variables[dep]

            if debug:
                print ' '
                print 'var =', var
                print 'varDepsDict.keys = ', varDepsDict.keys()
                print 'varDepsDict = ', varDepsDict
            # calculate the monitor variable
            value = self.calculate_MonitorVariable(var, varDepsDict)

            # Load the data into monitor variables
            varKind = monitorDefinition[var][0]
            if varKind == 'S':
                var_obj[n_step] = value
            if varKind == 'P':
                var_obj[n_step,:] = value
            if varKind == 'SL':
                var_obj[n_step,:] = value
            if varKind == '2D':
                var_obj[n_step,:] = value



        # Finished writing to monitor file.  Close it.
        monitor_file.close()

        # Close plasma_state
        plasma_state.close()

        print 'update_monitor_file finished'
        return 0


#----------------------------------------------------------------------------------------------
#
# calculate_MonitorVariable function
#
# For given montor variable it calculates the variable value from the netcdf variable objects
# in the varDepsDict dictionary on which the monitor variable depends.
#
#----------------------------------------------------------------------------------------------

    def calculate_MonitorVariable(self, var, varDepsDict):

        if debug:
            print 'calculate_MonitorVariable, var = ', var
            for dep_name in varDepsDict.keys():
                print dep_name, '= ', varDepsDict[dep_name][:]
                print dep_name, ' shape = ', varDepsDict[dep_name].shape
    
        # ------------- Thermal plasma species ________________
        if var == 'Te(0)':
            value = varDepsDict['Ts'][:][0][0]
    
        if var == 'Te_ave':
            Te = varDepsDict['Ts'][:][0]
            rho = varDepsDict['rho'][:]
            vol = varDepsDict['vol'][:]
            rho_eq = varDepsDict['rho_eq'][:]
            value = self.vol_average(Te, rho, vol, rho_eq)
    
        if var == 'Te(rho)':
            value = varDepsDict['Ts'][:][0]
    
        if var == 'Ti(0)':
            value = varDepsDict['Ti'][:][0]
    
        if var == 'Ti_ave':
            Ti = varDepsDict['Ti'][:]
            rho = varDepsDict['rho'][:]
            vol = varDepsDict['vol'][:]
            rho_eq = varDepsDict['rho_eq'][:]
            value = self.vol_average(Ti, rho, vol, rho_eq)
    
        if var == 'Ti(rho)':
            value = varDepsDict['Ti'][:]
    
        if var == 'ne(0)':
            value = varDepsDict['ns'][:][0][0]
    
        if var == 'ne_ave':
            ne = varDepsDict['ns'][:][0]
            rho = varDepsDict['rho'][:]
            vol = varDepsDict['vol'][:]
            rho_eq = varDepsDict['rho_eq'][:]
            value = self.vol_average(ne, rho, vol, rho_eq)
    
        if var == 'ne(rho)':
            value = varDepsDict['ns'][:][0]
    
        if var == 'n_H(0)':
            index = find_label_index_in_list(varDepsDict, 'S_name', 'H')
            value = varDepsDict['ns'][:][index][0]
    
        if var == 'n_H_ave':
            index = self.find_label_index_in_list(varDepsDict, 'S_name', 'H')
            n_H = varDepsDict['ns'][:][index]
            rho = varDepsDict['rho'][:]
            vol = varDepsDict['vol'][:]
            rho_eq = varDepsDict['rho_eq'][:]
            value = self.vol_average(n_H, rho, vol, rho_eq)
    
        if var == 'n_H(rho)':
            index = self.find_label_index_in_list(varDepsDict, 'S_name', 'H')
            value = varDepsDict['ns'][:][index]
    
        if var == 'n_D(0)':
            index = self.find_label_index_in_list(varDepsDict, 'S_name', 'D')
            value = varDepsDict['ns'][:][index][0]
    
        if var == 'n_D_ave':
            index = self.find_label_index_in_list(varDepsDict, 'S_name', 'D')
            n_D = varDepsDict['ns'][:][index]
            rho = varDepsDict['rho'][:]
            vol = varDepsDict['vol'][:]
            rho_eq = varDepsDict['rho_eq'][:]
            value = self.vol_average(n_D, rho, vol, rho_eq)
    
        if var == 'n_D(rho)':
            index = self.find_label_index_in_list(varDepsDict, 'S_name', 'D')
            value = varDepsDict['ns'][:][index]
    
        if var == 'n_T(0)':
            index = self.find_label_index_in_list(varDepsDict, 'S_name', 'T')
            value = varDepsDict['ns'][:][index][0]
    
        if var == 'n_T_ave':
            index = self.find_label_index_in_list(varDepsDict, 'S_name', 'T')
            n_T = varDepsDict['ns'][:][index]
            rho = varDepsDict['rho'][:]
            vol = varDepsDict['vol'][:]
            rho_eq = varDepsDict['rho_eq'][:]
            value = self.vol_average(n_T, rho, vol, rho_eq)
    
        if var == 'n_T(rho)':
            index = self.find_label_index_in_list(varDepsDict, 'S_name', 'T')
            value = varDepsDict['ns'][:][index]
    
        if var == 'n_He4(0)':
            index = self.find_label_index_in_list(varDepsDict, 'S_name', 'He4')
            value = varDepsDict['ns'][:][index][0]
    
        if var == 'n_He4_ave':
            index = self.find_label_index_in_list(varDepsDict, 'S_name', 'He4')
            n_He4 = varDepsDict['ns'][:][index]
            rho = varDepsDict['rho'][:]
            vol = varDepsDict['vol'][:]
            rho_eq = varDepsDict['rho_eq'][:]
            value = self.vol_average(n_He4, rho, vol, rho_eq)
    
        if var == 'n_He4(rho)':
            index = self.find_label_index_in_list(varDepsDict, 'S_name', 'He4')
            value = varDepsDict['ns'][:][index]
    
        if var == 'J_BS(rho)':
            curr_bootstrap = varDepsDict['curr_bootstrap'][:]
            rho = varDepsDict['rho'][:]
            area = varDepsDict['area'][:]
            rho_eq = varDepsDict['rho_eq'][:]
            value = self.zone_to_surf_dens(curr_bootstrap, rho, area, rho_eq)
    
        if var == 'I_BS(rho)':
            value = varDepsDict['curr_bootstrap'][:]
    
        if var == 'I_BS_cum(rho)':
            curr_bootstrap = varDepsDict['curr_bootstrap'][:]
            value = self.cumulate(curr_bootstrap)
    
        if var == 'I_BS_total':
            current = varDepsDict['curr_bootstrap'][:]
            value = sum(current)
    
        if var == 'J_OH(rho)':
            curr_ohmic = varDepsDict['curr_ohmic'][:]
            rho = varDepsDict['rho'][:]
            area = varDepsDict['area'][:]
            rho_eq = varDepsDict['rho_eq'][:]
            value = self.zone_to_surf_dens(curr_ohmic, rho, area, rho_eq)
    
        if var == 'I_OH(rho)':
            value = varDepsDict['curr_ohmic'][:]
    
        if var == 'I_OH_cum(rho)':
            current = varDepsDict['curr_ohmic'][:]
            value = self.cumulate(current)
    
        if var == 'I_OH_total':
            current = varDepsDict['curr_ohmic'][:]
            value = sum(current)
    
        if var == 'Pe_OH_dens(rho)':
            power = varDepsDict['pohme'][:]
            rho = varDepsDict['rho'][:]
            vol = varDepsDict['vol'][:]
            rho_eq = varDepsDict['rho_eq'][:]
            value = self.zone_to_surf_dens(power, rho, vol, rho_eq)
    
        if var == 'Pe_OH(rho)':
            value = varDepsDict['pohme'][:]

        if var == 'Pe_OH_cum(rho)':
            power = varDepsDict['pohme'][:]
            value = self.cumulate(power)
    
        if var == 'Pe_OH_total':
            power = varDepsDict['pohme'][:]
            value = sum(power)
    
        if var == 'Zeff(rho)':
            value = varDepsDict['Zeff'][:]
    
        if var == 'Zeff(0)':
            value = varDepsDict['Zeff'][:][0]
    
        if var == 'Zeff_ave':
            Zeff = varDepsDict['Zeff'][:]
            rho = varDepsDict['rho'][:]
            vol = varDepsDict['vol'][:]
            rho_eq = varDepsDict['rho_eq'][:]
            value = self.vol_average(Zeff, rho, vol, rho_eq)
    
    
    
        # ------------- Minority plasma species ________________
        if var == 'Pmin_e(rho)':
            value = varDepsDict['pmine'][:]
    
        if var == 'Pmin_i(rho)':
            value = varDepsDict['pmini'][:]
    
        if var == 'Eperp_mini(rho)':
            value = varDepsDict['eperp_mini'][:][0] # N.B. assumes one minority species
    
        if var == 'Epll_mini(rho)':
            value = varDepsDict['epll_mini'][:][0] # N.B. assumes one minority species
    
        # ------------- Equilibrium  ________________
    
        if var == 'P_eq(rho)':
            value = varDepsDict['P_eq'][:]
    
        if var == 'q_eq(rho)':
            value = varDepsDict['q_eq'][:]
    
        if var == 'q(0)':
            value = varDepsDict['q_eq'][:][0]
            # Large values of q make it hard to see more normal values, so suppress them.
            if value > 10.0: value = 10.0  
        if var == 'q(95)':
            q_eq = varDepsDict['q_eq'][:]
            rho_eq = varDepsDict['rho_eq'][:]
            psipol = varDepsDict['psipol'][:]
            value = self.calculate_qPsi(0.95, q_eq, rho_eq, psipol)
            # Large values of q make it hard to see more normal values, so suppress them.
            if value > 10.0: value = 10.0  
    
        if var == 'li(3)':
            curt = varDepsDict['curt'][:]
            q_eq = varDepsDict['q_eq'][:]
            grho2r2i = varDepsDict['grho2r2i'][:]
            phit = varDepsDict['phit'][:]
            R_axis = varDepsDict['R_axis'][:]
            vol = varDepsDict['vol'][:]
            value = self.calculate_li_3(curt, q_eq, grho2r2i, phit, R_axis, vol)

        if var == 'J_plasma(rho)':
            rho_eq = varDepsDict['rho_eq'][:]
            curt = varDepsDict['curt'][:]
            area = varDepsDict['area'][:]
            # curt is a cumulative variable, calcuate curt per zone
            curt_zone = []
            for i in range(1, len(curt)):
                curt_zone.append(curt[i] - curt[i-1])
            value = self.zone_to_surf_dens(curt_zone, rho_eq, area, rho_eq)  
    
        if var == 'I_plasma(rho)':
            curt = varDepsDict['curt'][:]
            value = []
            for i in range(len(curt)-1):
                value.append(curt[i+1] - curt[i])
            
        if var == 'I_plasma_cum(rho)':
            curt = varDepsDict['curt'][:]
            value = curt
    
        if var == 'I_plasma':
            value = varDepsDict['curt'][:][-1]
            
        if var == 'R_axis':
            value = varDepsDict['R_axis'][:]
            
        if var == 'Z_axis':
            value = varDepsDict['Z_axis'][:]
    
        if var == 'triang(95)':
            psipol = varDepsDict['psipol'][:]
            rho_eq = varDepsDict['rho_eq'][:]
            triang = varDepsDict['triang'][:]
            rho_eq_geo = varDepsDict['rho_eq_geo'][:]
            value = self.calculate_triangPsi(0.95, triang, psipol, rho_eq, rho_eq_geo)
    
        if var == 'elong(95)':
            psipol = varDepsDict['psipol'][:]
            rho_eq = varDepsDict['rho_eq'][:]
            elong = varDepsDict['elong'][:]
            rho_eq_geo = varDepsDict['rho_eq_geo'][:]
            value = self.calculate_elongPsi(0.95, elong, psipol, rho_eq, rho_eq_geo)
    
        if var == 'Psi(RZ)':
            value = varDepsDict['PsiRZ'][:]
    
        if var == 'rlim':
            value = varDepsDict['rlim'][:]
    
        if var == 'zlim':
            value = varDepsDict['zlim'][:]
            
        if var == 'N_GW':    
            ns = varDepsDict['ns'][:]
            rho = varDepsDict['rho'][:]
            rho_eq_geo = varDepsDict['rho_eq_geo'][:]
            R_surfMin = varDepsDict['R_surfMin'][:]
            R_surfMax = varDepsDict['R_surfMax'][:]
            curt = varDepsDict['curt'][:]
            value = self.calculate_N_GW(ns, rho, rho_eq_geo, R_surfMin, R_surfMax, curt)
    
        if var == 'beta_th':    
            P_eq = varDepsDict['P_eq'][:]
            rho = varDepsDict['rho'][:]
            rho_eq = varDepsDict['rho_eq'][:]
            vol = varDepsDict['vol'][:]
            B_axis_vac = varDepsDict['B_axis_vac'][:]
            mu0 = 4.0*math.pi/10000000.0
            value = mu0*self.vol_average(P_eq, rho, vol, rho_eq)/B_axis_vac**2
    
        if var == 'beta_N':    
            P_eq = varDepsDict['P_eq'][:]
            rho = varDepsDict['rho'][:]
            rho_eq = varDepsDict['rho_eq'][:]
            vol = varDepsDict['vol'][:]
            B_axis_vac = varDepsDict['B_axis_vac'][:]
            R_surfMin = varDepsDict['R_surfMin'][:]
            R_surfMax = varDepsDict['R_surfMax'][:]
            curt = varDepsDict['curt'][:]
            mu0 = 4.0*math.pi/10000000.0
            beta0 = mu0*self.vol_average(P_eq, rho, vol, rho_eq)/B_axis_vac**2
            a = (R_surfMax[-1] - R_surfMin[-1])/2.0
            value = 1.0e6*a*B_axis_vac*beta0/curt[-1]
    
        # ------------- Other State data  ________________
    
        if var == 'Vsurf':
            value = varDepsDict['vsur'][:]
    
        if var == 'V_loop(rho)':
            value = varDepsDict['V_loop'][:]
    
    
        # ------------- ECRF power deposition ________________
    
        if var == 'power_EC':
            power = varDepsDict['power_ec'][:]
            value = sum(power)
    
        if var == 'Pe_ecrf_dens(rho)':
            power = varDepsDict['peech'][:]
            rho = varDepsDict['rho_ecrf'][:]
            vol = varDepsDict['vol'][:]
            rho_eq = varDepsDict['rho_eq'][:]
            value = self.zone_to_surf_dens(power, rho, vol, rho_eq)
    
        if var == 'Pe_ecrf(rho)':
            value = varDepsDict['peech'][:]
 
        if var == 'Pe_ecrf_cum(rho)':
            power = varDepsDict['peech'][:]
            value = self.cumulate(power)

        if var == 'Pe_ecrf_total':
            power = varDepsDict['peech'][:]
            value = sum(power)
    
        if var == 'J_ecrf(rho)':
            curech = varDepsDict['curech'][:]
            rho = varDepsDict['rho_ecrf'][:]
            area = varDepsDict['area'][:]
            rho_eq = varDepsDict['rho_eq'][:]
            value = self.zone_to_surf_dens(curech, rho, area, rho_eq)
    
        if var == 'I_ecrf(rho)':
            value = varDepsDict['curech'][:]
    
        if var == 'I_ecrf_cum(rho)':
            curech = varDepsDict['curech'][:]
            value = self.cumulate(curech)
    
        if var == 'I_ecrf_total':
            current = varDepsDict['curech'][:]
            value = sum(current)
    
        # ------------- Lower Hybrid power deposition ________________
    
        if var == 'power_LH':
            power = varDepsDict['power_lh'][:]
            value = sum(power)
    
        if var == 'Pe_LH_dens(rho)':
            power = varDepsDict['pelh'][:]
            rho = varDepsDict['rho_lhrf'][:]
            vol = varDepsDict['vol'][:]
            rho_eq = varDepsDict['rho_eq'][:]
            value = self.zone_to_surf_dens(power, rho, vol, rho_eq)
    
        if var == 'Pe_LH(rho)':
            value = varDepsDict['pelh'][:]

        if var == 'Pe_LH_cum(rho)':
            power = varDepsDict['pelh'][:]
            value = self.cumulate(power)
    
        if var == 'Pe_LH_total':
            power = varDepsDict['pelh'][:]
            value = sum(power)
    
        if var == 'Pi_LH_dens(rho)':
            power = varDepsDict['pilh'][:]
            rho = varDepsDict['rho_lhrf'][:]
            vol = varDepsDict['vol'][:]
            rho_eq = varDepsDict['rho_eq'][:]
            value = self.zone_to_surf_dens(power, rho, vol, rho_eq)
    
        if var == 'Pi_LH(rho)':
            value = varDepsDict['pilh'][:]

        if var == 'Pi_LH_cum(rho)':
            power = varDepsDict['pilh'][:]
            value = self.cumulate(power)
    
        if var == 'Pi_LH_total':
            power = varDepsDict['pilh'][:]
            value = sum(power)
    
        if var == 'J_LH(rho)':
            curlh = varDepsDict['curlh'][:]
            rho = varDepsDict['rho_lhrf'][:]
            area = varDepsDict['area'][:]
            rho_eq = varDepsDict['rho_eq'][:]
            value = self.zone_to_surf_dens(curlh, rho, area, rho_eq)
    
        if var == 'I_LH(rho)':
            value = varDepsDict['curlh'][:]
    
        if var == 'I_LH_cum(rho)':
            curlh = varDepsDict['curlh'][:]
            value = self.cumulate(curlh)
    
        if var == 'I_LH_total':
            current = varDepsDict['curlh'][:]
            value = sum(current)
    
        # ------------- Fusion Product Fast Ions ________________
    
        if var == 'ni_FUS(rho)':
            ni_fus = varDepsDict['nfusi'][:]
            value = [0 for i in range(len(ni_fus[0])) ]
            for j in range(len(ni_fus)):
                for i in range(len(ni_fus[0])):
                    value[i] = value[i] + ni_fus[j][i]
    
        if var == 'Pe_FUS_dens(rho)':
            power = varDepsDict['pfuse'][:]
            rho = varDepsDict['rho_fus'][:]
            vol = varDepsDict['vol'][:]
            rho_eq = varDepsDict['rho_eq'][:]
            value = self.zone_to_surf_dens(power, rho, vol, rho_eq)
    
        if var == 'Pe_FUS(rho)':
            value = varDepsDict['pfuse'][:]

        if var == 'Pe_FUS_cum(rho)':
            power = varDepsDict['pfuse'][:]
            value = self.cumulate(power)
    
        if var == 'Pe_FUS_total':
            power = varDepsDict['pfuse'][:]
            value = sum(power)
    
        if var == 'Pi_FUS_dens(rho)':
            power = varDepsDict['pfusi'][:]
            rho = varDepsDict['rho_fus'][:]
            vol = varDepsDict['vol'][:]
            rho_eq = varDepsDict['rho_eq'][:]
            value = self.zone_to_surf_dens(power, rho, vol, rho_eq)
    
        if var == 'Pi_FUS(rho)':
            value = varDepsDict['pfusi'][:]

        if var == 'Pi_FUS_cum(rho)':
            power = varDepsDict['pfusi'][:]
            value = self.cumulate(power)
    
        if var == 'Pi_FUS_total':
            power = varDepsDict['pfusi'][:]
            value = sum(power)
    
        if var == 'J_FUS(rho)':
            curfusn = varDepsDict['curfusn'][:]
            rho = varDepsDict['rho_fus'][:]
            area = varDepsDict['area'][:]
            rho_eq = varDepsDict['rho_eq'][:]
            value = self.zone_to_surf_dens(curfusn, rho, area, rho_eq)
    
        if var == 'I_FUS(rho)':
            value = varDepsDict['curfusn'][:]
    
        if var == 'I_FUS_cum(rho)':
            curfusn = varDepsDict['curfusn'][:]
            value = self.cumulate(curfusn)
    
        if var == 'I_FUS_total':
            current = varDepsDict['curfusn'][:]
            value = sum(current)
    
        if var == 'Eperp_FUSI(rho)':
            eperp = varDepsDict['eperp_fusi'][:]
            value = [0 for i in range(len(eperp[0])) ]
            for j in range(len(eperp)):
                for i in range(len(eperp[0])):
                    value[i] = value[i] + eperp[j][i]
    
        if var == 'Epll_FUSI(rho)':
            epll = varDepsDict['epll_fusi'][:]
            value = [0 for i in range(len(epll[0])) ]
            for j in range(len(epll)):
                for i in range(len(epll[0])):
                    value[i] = value[i] + epll[j][i]
   
        if var == 'P_nuc_FUS':
            # A couple of things to note:
            #    This calls calculate_MonitorVariable( ) recursively
            #    This takes the instantaneous fusion power production to be 5 times
            #    the instantaneous plasma heating by fusion products.  It might not
            #    be such a good approximation for a machine smaller than ITER.
            Pe_FUS_total = self.calculate_MonitorVariable('Pe_FUS_total', varDepsDict)
            Pi_FUS_total = self.calculate_MonitorVariable('Pi_FUS_total', varDepsDict)
            value = 17.6/3.5*(Pe_FUS_total + Pi_FUS_total)
   
        if var == 'Q_nuc_FUS':
            # This calls calculate_MonitorVariable( ) recursively.
            
            P_nuc_FUS = self.calculate_MonitorVariable('P_nuc_FUS', varDepsDict)
            heating = self.calculate_MonitorVariable('Pe_OH_total', varDepsDict)
            
            if 'power_ec' in varDepsDict.keys():
	            heating = heating + self.calculate_MonitorVariable('power_EC', varDepsDict)
            if 'power_lh' in varDepsDict.keys():
	            heating = heating + self.calculate_MonitorVariable('power_LH', varDepsDict)
            if 'power_ic' in varDepsDict.keys():
	            heating = heating + self.calculate_MonitorVariable('power_IC', varDepsDict)
            if 'power_nbi' in varDepsDict.keys():
	            heating = heating + self.calculate_MonitorVariable('power_NB', varDepsDict)
            value = P_nuc_FUS/heating
    
        # ------------- ICRF power deposition ________________
    
        if var == 'power_IC':
            power = varDepsDict['power_ic'][:]
            value = sum(power)
    
        if var == 'Pe_icrf_dens(rho)':
            power = varDepsDict['picrf_totals'][:][0]
            rho = varDepsDict['rho_icrf'][:]
            vol = varDepsDict['vol'][:]
            rho_eq = varDepsDict['rho_eq'][:]
            value = self.zone_to_surf_dens(power, rho, vol, rho_eq)

        if var == 'Pe_icrf_cum(rho)':
            power = varDepsDict['picrf_totals'][:][0]
            value = self.cumulate(power)
    
        if var == 'Pe_icrf(rho)':
            value = varDepsDict['picrf_totals'][:][0]
    
        if var == 'Pe_icrf_total':
            power = varDepsDict['picrf_totals'][:][0]
            value = sum(power)
    
        if var == 'Pi_icrf_dens(rho)':
            power = varDepsDict['picth'][:]
            rho = varDepsDict['rho_icrf'][:]
            vol = varDepsDict['vol'][:]
            rho_eq = varDepsDict['rho_eq'][:]
            value = self.zone_to_surf_dens(power, rho, vol, rho_eq)
    
        if var == 'Pi_icrf(rho)':
            value = varDepsDict['picth'][:]

        if var == 'Pi_icrf_cum(rho)':
            power = varDepsDict['picth'][:]
            value = self.cumulate(power)
    
        if var == 'Pi_icrf_total':
            power = varDepsDict['picth'][:]
            value = sum(power)
    
        if var == 'J_icrf(rho)':
            curich = varDepsDict['curich'][:]
            rho = varDepsDict['rho_icrf'][:]
            area = varDepsDict['area'][:]
            rho_eq = varDepsDict['rho_eq'][:]
            value = self.zone_to_surf_dens(curich, rho, area, rho_eq)
   
        if var == 'I_icrf(rho)':
            value = varDepsDict['curich'][:]
    
        if var == 'I_icrf_cum(rho)':
            curich = varDepsDict['curich'][:]
            value = self.cumulate(curich)
     
        if var == 'I_icrf_total':
            current = varDepsDict['curich'][:]
            value = sum(current)
    
        # ------------- Minority plasma species ________________
    
        if var == 'nmin_icrf(rho)':
            ni_min = varDepsDict['nmini'][:]
            value = [0 for i in range(len(ni_min[0])) ]
            for j in range(len(ni_min)):
                for i in range(len(ni_min[0])):
                    value[i] = value[i] + ni_min[j][i]
    
        if var == 'Pmin_e_dens(rho)':
            power = varDepsDict['pmine'][:]
            rho = varDepsDict['rho_icrf'][:]
            vol = varDepsDict['vol'][:]
            rho_eq = varDepsDict['rho_eq'][:]
            value = self.zone_to_surf_dens(power, rho, vol, rho_eq)
    
        if var == 'Pmin_e(rho)':
            value = varDepsDict['pmine'][:]
 
        if var == 'Pmin_e_cum(rho)':
            power = varDepsDict['pmine'][:]
            value = self.cumulate(power)
    
        if var == 'Pmin_e_total':
            power = varDepsDict['pmine'][:]
            value = sum(power)

        if var == 'Pmin_i_dens(rho)':
            power = varDepsDict['pmini'][:]
            rho = varDepsDict['rho_icrf'][:]
            vol = varDepsDict['vol'][:]
            rho_eq = varDepsDict['rho_eq'][:]
            value = self.zone_to_surf_dens(power, rho, vol, rho_eq)
    
        if var == 'Pmin_i(rho)':
            value = varDepsDict['pmini'][:]
 
        if var == 'Pmin_i(rho)_cum':
            power = varDepsDict['pmini'][:]
            value = self.cumulate(power)
    
        if var == 'Pmin_i_total':
            power = varDepsDict['pmini'][:]
            value = sum(power)
    
        if var == 'Eperp_mini(rho)':
            eperp = varDepsDict['eperp_mini'][:]
            value = [0 for i in range(len(eperp[0])) ]
            for j in range(len(eperp)):
                for i in range(len(eperp[0])):
                    value[i] = value[i] + eperp[j][i]
    
        if var == 'Epll_mini(rho)':
            epll = varDepsDict['epll_mini'][:]
            value = [0 for i in range(len(epll[0])) ]
            for j in range(len(epll)):
                for i in range(len(epll[0])):
                    value[i] = value[i] + epll[j][i]
    
        # ------------- NBI power deposition ________________
    
        if var == 'power_NB':
            power = varDepsDict['power_nbi'][:]
            value = sum(power)
    
        # ------------- Neutral Beam plasma species ________________
    
        if var == 'ni_NB(rho)':
            ni_nb = varDepsDict['nbeami'][:]
            value = [0 for i in range(len(ni_nb[0])) ]
            for j in range(len(ni_nb)):
                for i in range(len(ni_nb[0])):
                    value[i] = value[i] + ni_nb[j][i]
    
        if var == 'Pe_NB_dens(rho)':
            power = varDepsDict['pbe'][:]
            rho = varDepsDict['rho_nbi'][:]
            vol = varDepsDict['vol'][:]
            rho_eq = varDepsDict['rho_eq'][:]
            value = self.zone_to_surf_dens(power, rho, vol, rho_eq)
    
        if var == 'Pe_NB(rho)':
            value = varDepsDict['pbe'][:]

        if var == 'Pe_NB_cum(rho)':
            power = varDepsDict['pbe'][:]
            value = self.cumulate(power)
    
        if var == 'Pe_NB_total':
            power = varDepsDict['pbe'][:]
            value = sum(power)
    
        if var == 'Pi_NB_dens(rho)':
            power = varDepsDict['pbi'][:]
            rho = varDepsDict['rho_nbi'][:]
            vol = varDepsDict['vol'][:]
            rho_eq = varDepsDict['rho_eq'][:]
            value = self.zone_to_surf_dens(power, rho, vol, rho_eq)
    
        if var == 'Pi_NB(rho)':
            value = varDepsDict['pbi'][:]

        if var == 'Pi_NB_cum(rho)':
            power = varDepsDict['pbi'][:]
            value = self.cumulate(power)
    
        if var == 'Pi_NB_total':
            power = varDepsDict['pbi'][:]
            value = sum(power)
    
        if var == 'J_NB(rho)':
            curbeam = varDepsDict['curbeam'][:]
            rho = varDepsDict['rho_nbi'][:]
            area = varDepsDict['area'][:]
            rho_eq = varDepsDict['rho_eq'][:]
            value = self.zone_to_surf_dens(curbeam, rho, area, rho_eq)
    
        if var == 'I_NB(rho)':
            value = varDepsDict['curbeam'][:]
    
        if var == 'I_NB_cum(rho)':
            curbeam = varDepsDict['curbeam'][:]
            value = self.cumulate(curbeam)
    
        if var == 'I_NB_total':
            current = varDepsDict['curbeam'][:]
            value = sum(current)
    
        if var == 'Eperp_NBI(rho)':
            eperp = varDepsDict['eperp_beami'][:]
            value = [0 for i in range(len(eperp[0])) ]
            for j in range(len(eperp)):
                for i in range(len(eperp[0])):
                    value[i] = value[i] + eperp[j][i]
    
        if var == 'Epll_NBI(rho)':
            epll = varDepsDict['epll_beami'][:]
            value = [0 for i in range(len(epll[0])) ]
            for j in range(len(epll)):
                for i in range(len(epll[0])):
                    value[i] = value[i] + epll[j][i]
    
        if debug:
            print ' '
            print var, ' = ', value
    
        return value
        
        
# ------------------------------------------------------------------------------
#
# "Private"  methods
#
# ------------------------------------------------------------------------------

    # Try to get config parameter - wraps the exception handling for get_config_parameter()
    def get_config_param(self, services, param_name, optional=False):

        try:
            value = services.get_config_param(param_name)
            print param_name, ' = ', value
        except Exception :
            if optional: 
                print 'config parameter ', param_name, ' not found'
                value = None
            else:
                message = 'required config parameter ', param_name, ' not found'
                print message
                services.exception(message)
                raise
        
        return value

    # Try to get component specific config parameter - wraps the exception handling
    def get_component_param(self, services, param_name, optional=False):

        if hasattr(self, param_name):
            value = getattr(self, param_name)
            print param_name, ' = ', value
        elif optional:
            print 'optional config parameter ', param_name, ' not found'
            value = None
        else:
            message = 'required component config parameter ', param_name, ' not found'
            print message
            services.exception(message)
            raise
        
        return value

