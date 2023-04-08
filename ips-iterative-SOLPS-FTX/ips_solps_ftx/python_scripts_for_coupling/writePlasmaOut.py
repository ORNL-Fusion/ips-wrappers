#!/usr/bin/env python
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import pickle
from fort44ext import fort44ext # fort.44 parser
from b2fextract import b2fextract # b2fstate parser
from rizp_extract import rizp_extract # return rizp from ionization_potentials.txt

# The purpose of this code is to calculate the target heat flux from SOLPS-ITER without relying on the b2plot wlld script, but instead using the standard output of SOLPS-ITER.
# read and post-process plasma heat load from b2fstate
# read and post-process neutral heat load from fort.44


### Example of input parameter for the funcion ###
#rad_grid = [17,19,20,22,25] # index convention? Fortran (Python) or MATLAB?

def SOLPS_heatflux_for_FTX(rad_grid=1, b2fstate_file='b2fstate', b2fgmtry_file='b2fgmtry', fort44_file='fort.44', outputFile=None, print_test=False, logFile=None): 

    cwd = os.getcwd()
    print(' ')
    print('\t Called SOLPS_heatflux_for_FTX:')

    #if pikle file exists, read from pkl file:
    #more elegant: turn into loop over keys,values
    pkl_file=cwd+'/writePlasmaOut.pkl'
    if os.path.exists(pkl_file):
        print('\t use inputs from pkl file: ', pkl_file)
        dic = pickle.load( open( pkl_file, "rb" ) )
        #first check the log file, to print everything there
        if 'logFile' in dic:
            logFile=dic['logFile']
        rad_grid=dic['rad_grid']
        b2fstate_file=dic['b2fstate_file']
        b2fgmtry_file=dic['b2fgmtry_file']
        fort44_file=dic['fort44_file']
        print_test=dic['print_test']        
        outputFile=dic['outputFile']
    else:
        print('no pkl file found, continue with function-call-inputs or defaults')

    if print_test:
        print(' ')
        print('\t launched script with inputs:')
        print('\t rad_grid = ', rad_grid)
        print('\t b2fstate_file = ', b2fstate_file)
        print('\t b2fgmtry_file = ', b2fgmtry_file)
        print('\t fort44_file = ', fort44_file)
        print('\t logFile =', logFile)        

    print(' ')
    if logFile is not None:
        print('\t printing to log file defined in keywords: ')
        print('\t', logFile)
        sys.stdout.flush()
        logF=open(logFile , 'a')
        orig_sys = sys.stdout
        sys.stdout = logF
    else:
        print('\t No log file defined; using default sys.stdout')
        sys.stdout.flush()
        
    if outputFile is not None:
        print('\t writing output to file defined in keywords: ')
        print('\t', outputFile)
        outFile=open(outputFile , 'a')
    elif logFile is not None:
        print('\t No output file defined, but log file exists. Write output to:')
        print('\t', logFile)
        outFile=logF
    else:
        print('\t No log or output files defined; using default sys.stdout')
    sys.stdout.flush()
    
    # constant
    ec = 1.60217662e-19
    
    # Read grid dimension and species info from fort.44
    with open(fort44_file, "r") as f:
        # Read the first line
        first_line = f.readline()
        values = first_line.split()
        nx = int(values[0])
        ny = int(values[1])
        pol_B25 = nx + 2
        rad_B25 = ny + 2
        # Read the second line
        second_line = f.readline()
        values = second_line.split()
        natom = int(values[0])
        nmol = int(values[1])
        nion = int(values[2])
    
        # Read and process species names
        species_names = []
        for _ in range(natom):
            species_line = f.readline()
            species = species_line.strip()
            formatted_species = species[0].upper() + species[1:].lower()
            species_names.append(formatted_species)
    
        # Join species names with '+'
        species_string = '+'.join(species_names)
    print('species = ',species_string)
    
    # Load rizp corresponding to species_string from ionization_potentials.txt
    rizp = rizp_extract(species_string)  # gives ionization potential in (eV) with size ns
    
    # Load geometry variables from b2fgmtry
    dv = b2fextract('vol',b2fgmtry_file) # cell volume
    gs = b2fextract('gs',b2fgmtry_file) # area
    gsx = gs[:,:,0] # area, x direction (target area etc.)


    # Plasma heat load calculation    
    # Load variables from b2fstate
    zamax = b2fextract('zamax',b2fstate_file,'no_reshape')
    ns = len(zamax)
    te = b2fextract('te',b2fstate_file)/ec  # in [eV] (by /ec)
    ti = b2fextract('ti',b2fstate_file)/ec  # in [eV] (by /ec)
    kinrgy = b2fextract('kinrgy',b2fstate_file)/ec # in [eV] (by /ec)
    ne = b2fextract('ne',b2fstate_file)
    #ua = b2fextract('ua',b2fstate_file)
    na = b2fextract('na',b2fstate_file)
    po = b2fextract('po',b2fstate_file) # in [V]
    # there is no fne only _32 and _52
    #fne_32 = b2fextract('fne_32',b2fstate_file)
    #fne_52 = b2fextract('fne_52',b2fstate_file)
    #fne_32_52 = fne_32+fne_52
    fna = b2fextract('fna',b2fstate_file) # fna = fna_32+fna_52 (checked 20230317)
    #fna_32 = b2fextract('fna_32',b2fstate_file)
    #fna_52 = b2fextract('fna_52',b2fstate_file)
    fch = b2fextract('fch',b2fstate_file) # in [A]
    fhe = b2fextract('fhe',b2fstate_file)
    fhi = b2fextract('fhi',b2fstate_file)
    
    # Split fluxes to x- (pol) and y- (rad) direction to remove [x,y] dimension
    fhex = fhe[:,:,0]
    fhey = fhe[:,:,1]
    fhix = fhi[:,:,0]
    fhiy = fhi[:,:,1]
    #fnex_32_52 = fne_32_52[:,:,0]
    #fney_32_52 = fne_32_52[:,:,1]
    fnax = fna[:,:,::2]  # select 0,2,4,6,... (nx+2)*(ny+2)*ns
    fnay = fna[:,:,1::2]  # select 1,3,5,7,.. (nx+2)*(ny+2)*ns
    fchx = fch[:,:,0]
    temp = np.zeros_like(fchx)
    for is_ in range(ns):
        temp += zamax[is_] * fnax[:, :, is_]
    fnex = temp - fchx/ec # different from fnex_32+fnex_52      
    
    
    # map cell centered quantities to cell face using average (1/2 weight) as wlld does
    te_face = np.zeros((pol_B25, rad_B25))
    ti_face = np.zeros((pol_B25, rad_B25))
    po_face = np.zeros((pol_B25, rad_B25))
    kinrgy_face = np.zeros((pol_B25, rad_B25, ns))
    for iy in range(rad_B25):
        for ix in range(1, pol_B25):
            te_face[ix, iy] = 0.5 * (te[ix - 1, iy] + te[ix, iy])
            ti_face[ix, iy] = 0.5 * (ti[ix - 1, iy] + ti[ix, iy])
            for is_ in range(ns):
                kinrgy_face[ix, iy, is_] = 0.5 * (kinrgy[ix - 1, iy, is_] + kinrgy[ix, iy, is_])
            po_face[ix, iy] = 0.5 * (po[ix - 1, iy] + po[ix, iy])
    
    
    # sum over species for fna related quantities that has extra ns dimension
    fnax_energy = np.zeros((pol_B25, rad_B25))  # sum over species
    for is_ in range(ns):
        #fnax_energy += np.squeeze(fnax[:, :, is_]) * (ti_face + np.squeeze(kinrgy_face[:, :, is_]) + rizp[is_]) * ec
        fnax_energy += fnax[:, :, is_] * (ti_face + kinrgy_face[:, :, is_] + rizp[is_]) * ec
    
    
    # plasma contribution, only x direction
    total_energy_flux = fhex + fhix \
        + fnex * te_face * ec \
        + fnax_energy
 #       - fchx * po_face # [Ampere*Volt]
    
    
    
    
    # Neutral heat load calculation
    # Read _res info from fort.44
    eirdiag = fort44ext('eirdiag',fort44_file) # info linking _res indices to (ix,iy)
    # atom
    ewlda_res = fort44ext('ewlda_res',fort44_file)
    ewldea_res = fort44ext('ewldea_res',fort44_file)
    # mol
    ewldm_res = fort44ext('ewldm_res',fort44_file)
    ewldem_res = fort44ext('ewldem_res',fort44_file)
    ewldmr_res = fort44ext('ewldmr_res',fort44_file) # recombination of atom to mol (nmol*NCL)
    # etc.
    ewldrp_res = fort44ext('ewldrp_res',fort44_file) # kinrgy carried by recycling atom and mol (NCL)
    sarea_res = fort44ext('sarea_res',fort44_file) # area (_res quantities are already divided by area so W/m^2 unit)
    
    
    # reshape _res into natom (nmol) * nsurf and sum over natom (nmol)
    length_atom = len(ewlda_res)
    nsurf_atom = length_atom // natom
    
    length_mol = len(ewldm_res)
    nsurf_mol = length_mol // nmol
    
    ewlda_res_2d = ewlda_res.reshape(natom, nsurf_atom, order = 'F')
    summed_ewlda_res = ewlda_res_2d.sum(axis=0)
    
    ewldea_res_2d = ewldea_res.reshape(natom, nsurf_atom, order = 'F')
    summed_ewldea_res = ewldea_res_2d.sum(axis=0)
    
    ewldm_res_2d = ewldm_res.reshape(nmol, nsurf_mol, order = 'F')
    summed_ewldm_res = ewldm_res_2d.sum(axis=0)
    
    ewldem_res_2d = ewldem_res.reshape(nmol, nsurf_mol, order = 'F')
    summed_ewldem_res = ewldem_res_2d.sum(axis=0)
    
    ewldmr_res_2d = ewldmr_res.reshape(nmol, nsurf_mol, order = 'F')
    summed_ewldmr_res = ewldmr_res_2d.sum(axis=0)
    
    
    
    # total neutral heat load: sum over _res quantities (now size matched)
    ewld_load_sum = summed_ewlda_res + summed_ewldm_res - summed_ewldea_res - summed_ewldem_res + summed_ewldmr_res # well matched! ewldrp goes to plasma contribution
    ewld_load_sum = ewld_load_sum
    
    
    # Read indices for _res and corresponding (pol,rad)
    # Reshape eirdiag from 1 dim array to 2 dim
    # eirdiag is a 1-dimensional array of shape (5 * n + 1,) in fort.44
    n = len(eirdiag) - 1 # make it 5*n form
    num_rows = 5 # eirdiag fixed format, always 5 rows in fort.44
    first_row_size = n // num_rows + 1 # first row contains additional column with the #total additional suface
    other_row_size = n // num_rows # will add 0 to the other rows to make it full 2d array with size below
    eirdiag_reshaped = np.zeros((num_rows, first_row_size)) # initialize
    eirdiag_reshaped[0, 0:first_row_size] = eirdiag[0:first_row_size]
    start_idx = first_row_size
    for i in range(1, num_rows):
        eirdiag_reshaped[i, 0:other_row_size] = eirdiag[start_idx:start_idx + other_row_size]
        start_idx += other_row_size
    # Now eirdiag_reshaped is a num_rows x first_row_size array as desired
    for i in range(0,5):
        #print(eirdiag_reshaped[i,:])
    # Caution! each column gives each NSTS but it is not ordered in the same way as in the input.dat
    
    # identify eirdiag: find radial surface (e.g., target)
    # and find its (pol,rad) indices and corresponding _res array indices
        num_columns = eirdiag_reshaped.shape[1]
        res_indices = [] # to index _res parameters 1 dimensional
        pr_indices = [] # to index total_energy_flux (pol,rad) 2 dimensional
    
    for col in range(num_columns - 1):  # exclude last column since it is dummy (for additional surf)
        res_start = int(eirdiag_reshaped[0, col])
        res_end = int(eirdiag_reshaped[0, col + 1])
        res_indices.append((res_start, res_end - 1))
    
        surface_type = int(eirdiag_reshaped[1, col])
        surface_index = int(eirdiag_reshaped[2, col])
        surface_start = int(eirdiag_reshaped[3, col])
        surface_end = int(eirdiag_reshaped[4, col])
    
        if surface_type == 1:  # Poloidal surface
            pr_indices.append((slice(surface_start, surface_end + 1), surface_index)) # end+1 for python indexing rule
        elif surface_type == 2:  # Radial surface
            pr_indices.append((surface_index, slice(surface_start, surface_end + 1)))
    # Now res_indices contains the _res index range for each column,
    # pr_indices shows (poloidal, radial) index for each column.
    
    res_slices = [slice(start, end) for start, end in res_indices] # gives indexing info for _res
    # surface_start:surface_end+1, res_start:res_end+1 does not match, 1 dummy index at _res end -> sum or whatever
    
    
    
    neutral_energy_flux = np.zeros((pol_B25, rad_B25))
    reflected_energy_flux = np.zeros((pol_B25, rad_B25))
    for col, surface_type in enumerate(eirdiag_reshaped[1, :]):
        if surface_type == 2:  # Radial surface
            #print('pr_indices[col]',pr_indices[col])
            #print('res_slices[col]',res_slices[col])
            neutral_energy_flux[pr_indices[col][0]+1, pr_indices[col][1]] += ewld_load_sum[res_slices[col]]
            reflected_energy_flux[pr_indices[col][0]+1, pr_indices[col][1]] += ewldrp_res[res_slices[col]]
            # +1 to pol coordinate to match flux convention in B2.5 (always left cell surface)
    
   
    
    # Wpls: divide by area to make it [W/m^2], subtract reflected energy
    if print_test:
        print('check value of area (non zero) before deviding energy flux')
        print('gsx = ', gsx)
        sys.stdout.flush()
    total_energy_flux_outer = total_energy_flux/gsx - reflected_energy_flux
    total_energy_flux_inner = -total_energy_flux/gsx - reflected_energy_flux

    # Wpls + Wneut
    total_energy_flux_add_neutral_inner = total_energy_flux_inner + neutral_energy_flux
    total_energy_flux_add_neutral_outer = total_energy_flux_outer + neutral_energy_flux

    if print_test:
        print(' ')
        print('\t check where to return the output and close files')

    if (logFile is not None) or (outputFile is not None): 
        outFile.write(str(total_energy_flux_add_neutral_outer[pol_B25-1,rad_grid]))
        if print_test:
            print('\t wrote result into ', outFile.name)

        try:
            if print_test:
                print('\t close outputFile', outFile.name)
            outFile.close()
        except:
            if print_test:
                print('\t could not close outputFile')
                
        if logFile is not None:
            if print_test:
                print('\t close logFile ', logF.name)
            sys.stdout.flush()
            logF.close()
            sys.stdout = orig_sys

        else:
            if print_test:
                print('\t could not close logFile')
        return
    else:
        print('\t no log or output file. return output')
        sys.stdout.flush()
        return total_energy_flux_add_neutral_outer[pol_B25-1,rad_grid]
    
if __name__ == '__main__':

   import shutil

   SOLPS_heatflux_for_FTX()
