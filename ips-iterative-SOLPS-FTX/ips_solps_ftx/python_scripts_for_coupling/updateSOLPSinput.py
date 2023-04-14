#!/usr/bin/env python
import numpy as np
import os
import matplotlib.pyplot as plt
from b2fextract import b2fextract
from fort44ext import fort44ext # fort.44 parser

# Jae-Sun Park
# parkj@ornl.gov
# 2023

# 1) read Twall, 'R's from ftxOut.txt
# 2) READ TRIM database and suggest scale factor RECYCF
# 3) change EWALL, RECYCF, RECYCT in input.dat

# Requires manual setup of target surface elements in input.dat block 3a (SOLPS)
# Provide SURFMOD_names (rad_grid_name below in USER DEFINED INPUTS section) correspond to grid # (in ftxOut.txt)


### USER DEFINED INPUTS ###
# SOLPS run directory
cwd = os.getcwd()
directory = '/global/homes/p/parkjs/SOLPS-FTX_atom/coupling_IO/' # path for input files include '/' at the end
dirlist = ['run_new_uinp_NTCPU_20','run_new_uinp_NTCPU_50','run_new_uinp_NTCPU_150'] # if several subdirectory exists, use this
TRIM_file = directory+'D_on_W'  # TRIM file to calculate RF_TRIM
FTX_file = cwd+'/ftxOut.txt'  # FTX output file
inputdat_file = cwd+'/input.dat'
fort44_file = cwd+'/fort.44'
rad_grid_name = ['Left1', 'DIMES', 'Right1'] # corresponds to grid # in ftxOut.txt matches with SURFMOD_name in input.dat
b2fstate_file = cwd+'/b2fstate' # to get Te, Ti info
b2fgmtry_file = cwd+'/b2fgmtry' # to get B field info
ec = 1.60217662e-19 # electron charge in [C]
###########################

###################
### RUN SCRIPTS ###
### after defs ####
###################

def read_lines(filename):
    try:
        f = open(filename,'r')
    except:
        print("No %s file" % (filename))
    temp = f.readlines()
    f.close()
    return temp



# Fine line number containing given word
def find_line(lines,word):
    line_num = -1
    for line in lines:
        line_num += 1
        if word in line:
            return line_num

# replace string
def replace_str_index(text,index_init,index_end,replacement=''):
    return '%s%s%s'%(text[:index_init],replacement,text[index_end+1:])



def get_FTX_out(FTXfile):
    lines = read_lines(FTXfile)
    rad_grid = []
    Twall = []
    R_FT = []
    R_X = []
    R_FTX = []
    for i in range(1,len(lines)): # start from the 2nd line (data body)
        rad_grid.append(lines[i].split()[0])
        Twall.append(lines[i].split()[1])
        R_FT.append(lines[i].split()[2])
        R_X.append(lines[i].split()[3])
        R_FTX.append(lines[i].split()[4]) # = R_FT + (1-R_FT)*R_X
    return rad_grid, Twall, R_FT, R_X, R_FTX
    




###### RECYCF: new approach based on wld #####
#FTX_file = directory+'ftxOut.txt'  # FTX output file
#inputdat_file = directory+'input.dat'
#rad_grid_name = ['Left1', 'DIMES', 'Right1'] # corresponds to grid # in ftxOut.txt matches with SURFMOD_name in input.dat



def find_numstrata_in_block14(lines, start_line, end_line):
    for i in range(start_line, end_line):
        line = lines[i].strip()
        if len(line.split()) == 1 and line.split()[0].isdigit():
            arbitrary_number = int(line.split()[0])
            next_lines_numbers = []
            count = 0
            num_lines = ((arbitrary_number) // 12) + 1

            for j in range(1, num_lines + 1):
                numbers_in_line = [int(x) for x in lines[i + j].split()]
                next_lines_numbers.extend(numbers_in_line)
                count += len(numbers_in_line)
                print('numbers_in_line',numbers_in_line)
                print('next_lines_numbers',next_lines_numbers)

            if sum(next_lines_numbers) == arbitrary_number:
                return arbitrary_number, i + num_lines + 1

    return -1





"""
lines = read_lines(inputdat_file)
line_block_3a = find_line(lines,'Data for non-default standard surfaces')
line_block_3b = find_line(lines, 'Data for additional surfaces')
line_block_14 = find_line(lines,'Data for interfacing routine "infusr"')
line_block_15 = find_line(lines, 'Data for interfacing routine "geousr"')

nstra, end_line_numstrata_block = find_numstrata_in_block14(lines, line_block_14, line_block_15)

nlim = int(lines[line_block_3b+1].split()[0])
print('nlim',nlim)

if end_line_numstrata_block != -1:
    print(f"Numstrata block ends at line: {end_line_numstrata_block}")
else:
    print("Numstrata block in not found.")

# get info from fort.44
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


wldra = fort44ext('wldra(0)',fort44_file) # atom -> atom
wldpa = fort44ext('wldpa(0)',fort44_file) # bulk ion -> atom
wldna = fort44ext('wldna(0)',fort44_file) # atom incidence
wldpp = fort44ext('wldpp(0)',fort44_file) # bluk ion incidence, different dim with above 3 data
# dimension = (NLIM + NSTS + 1 (dummy), NATM (NPLS for wldpp))

length_atom = len(wldra)
nsurf_atom = length_atom // natom # NLIM + NSTS + 1

length_pls = len(wldpp)
npls = length_pls // nsurf_atom


wldra_2d = wldra.reshape(nsurf_atom, natom, order = 'F')
wldpa_2d = wldpa.reshape(nsurf_atom, natom, order = 'F')
wldna_2d = wldna.reshape(nsurf_atom, natom, order = 'F')
wldpp_2d = wldpp.reshape(nsurf_atom, npls, order = 'F')

#print('wldra',wldra_2d[66,0])
#print('wldpa',wldpa_2d[66,0])
#print('wldna',wldna_2d[66,0])
#print('wldpp',wldpp_2d[66,0])




for surfmod_name in rad_grid_name:
    print(surfmod_name)
    line_SURFMOD_def = find_line(lines[line_block_3a:line_block_3b], 'SURFMOD_'+surfmod_name)
    # former 3 line gives info of that strata (-2 contains geometry information)
#    print(lines[line_block_3a+line_SURFMOD_def]) # SURFMOD line
    pol_location = int(lines[line_block_3a+line_SURFMOD_def-2].split()[2])-1
    rad_location_1 = lines[line_block_3a+line_SURFMOD_def-2].split()[3]
    rad_location_2 = lines[line_block_3a+line_SURFMOD_def-2].split()[4]
    print('pol_location',pol_location)
    print('rad_location_1',rad_location_1)
    print('rad_location_2',rad_location_2)

    matching_line = -1

    print("end_line_numstrata_block + 1", end_line_numstrata_block + 1)
    print("end_line_numstrata_block + 1 + nstra", end_line_numstrata_block + 1 + nstra)

    for i in range(end_line_numstrata_block + 1, end_line_numstrata_block + 1 + nstra):
        elements = lines[i].split()
        if (elements[1] == str(pol_location)) and (elements[4] == str(rad_location_1)) and (elements[5] == str(rad_location_2)):
            matching_line = i - end_line_numstrata_block +1
            break

    if matching_line != -1:
        print("strata number: ", matching_line)
    else:
        print("cannnot find strata with given pol, rad info")

    R_F_SOLPS = (wldra_2d[nlim+matching_line-1,0]+wldpa_2d[nlim+matching_line-1,0])/(wldna_2d[nlim+matching_line-1,0]+wldpp_2d[nlim+matching_line-1,0])
    print('R_F_SOLPS',R_F_SOLPS)
"""
   



def calc_RECYCF(inputdat_file, fort44_file, rad_grid_name, R_FT):

    lines = read_lines(inputdat_file)
    line_block_3a = find_line(lines,'Data for non-default standard surfaces')
    line_block_3b = find_line(lines, 'Data for additional surfaces')
    line_block_14 = find_line(lines,'Data for interfacing routine "infusr"')
    line_block_15 = find_line(lines, 'Data for interfacing routine "geousr"')
    
    nstra, end_line_numstrata_block = find_numstrata_in_block14(lines, line_block_14, line_block_15)
    
    nlim = int(lines[line_block_3b+1].split()[0])
    print('nlim',nlim)
    
    if end_line_numstrata_block != -1:
        print(f"Numstrata block ends at line: {end_line_numstrata_block}")
    else:
        print("Numstrata block in not found.")
    
    # get info from fort.44
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
    
    
    wldra = fort44ext('wldra(0)',fort44_file) # atom -> atom
    wldpa = fort44ext('wldpa(0)',fort44_file) # bulk ion -> atom
    wldna = fort44ext('wldna(0)',fort44_file) # atom incidence
    wldpp = fort44ext('wldpp(0)',fort44_file) # bluk ion incidence, different dim with above 3 data
    # dimension = (NLIM + NSTS + 1 (dummy), NATM (NPLS for wldpp))
    
    length_atom = len(wldra)
    nsurf_atom = length_atom // natom # NLIM + NSTS + 1
    
    length_pls = len(wldpp)
    npls = length_pls // nsurf_atom
    
    
    wldra_2d = wldra.reshape(nsurf_atom, natom, order = 'F')
    wldpa_2d = wldpa.reshape(nsurf_atom, natom, order = 'F')
    wldna_2d = wldna.reshape(nsurf_atom, natom, order = 'F')
    wldpp_2d = wldpp.reshape(nsurf_atom, npls, order = 'F')
    
    #print('wldra',wldra_2d[66,0])
    #print('wldpa',wldpa_2d[66,0])
    #print('wldna',wldna_2d[66,0])
    #print('wldpp',wldpp_2d[66,0])         
    
    R_FT_counter = 0
    RECYCF = []
    for surfmod_name in rad_grid_name:
        print(surfmod_name)
        line_SURFMOD_def = find_line(lines[line_block_3a:line_block_3b], 'SURFMOD_'+surfmod_name)
        # former 3 line gives info of that strata (-2 contains geometry information)
    #    print(lines[line_block_3a+line_SURFMOD_def]) # SURFMOD line
        pol_location = int(lines[line_block_3a+line_SURFMOD_def-2].split()[2])-1
        rad_location_1 = lines[line_block_3a+line_SURFMOD_def-2].split()[3]
        rad_location_2 = lines[line_block_3a+line_SURFMOD_def-2].split()[4]
        print('pol_location',pol_location)
        print('rad_location_1',rad_location_1)
        print('rad_location_2',rad_location_2)
    
        matching_line = -1
    
        print("end_line_numstrata_block + 1", end_line_numstrata_block + 1)
        print("end_line_numstrata_block + 1 + nstra", end_line_numstrata_block + 1 + nstra)
    
    
        for i in range(end_line_numstrata_block + 1, end_line_numstrata_block + 1 + nstra):
            elements = lines[i].split()
            if (elements[1] == str(pol_location)) and (elements[4] == str(rad_location_1)) and (elements[5] == str(rad_location_2)):
                matching_line = i - end_line_numstrata_block +1
                break
    
        if matching_line != -1:
            print("strata number: ", matching_line)
        else:
            print("cannnot find strata with given pol, rad info")
    
        R_F_SOLPS = (wldra_2d[nlim+matching_line-1,0]+wldpa_2d[nlim+matching_line-1,0])/(wldna_2d[nlim+matching_line-1,0]+wldpp_2d[nlim+matching_line-1,0])
        print('R_F_SOLPS',R_F_SOLPS)
        R_FT_this = float(R_FT[R_FT_counter])
        print('R_FT',R_FT_this)

        RECYCF_this = R_FT_this/R_F_SOLPS
        print('RECYCF',RECYCF_this)

        RECYCF.append(RECYCF_this)

        R_FT_counter += 1

    return RECYCF
        













################################## <WORKING ###################################

# get Te and Ti at rad_OT from b2fstate
# use my own b2f parser b2fextract.py
# translate it to E0_mean


# 1) translate plasma information (Te, Ti, Bin) into (E0_mean, Ain_Mean)
#    input: b2fstate, b2fgmtry; output: E0_mean, Ain_mean
# 2) Generate samples of [E0, Ain] with [E0_mean, Ain_mean] maxwellian
# 3) 2D interpolation of each samples and get distributin of R_F_TRIM
# 4) Get averaged R_F_TRIM
# 5) RECYCF = R_FT/R_F_TRIM



def get_Te_Ti_Bin(rad_grid, b2fstate_file, b2fgmtry_file):
    # Assume outer target (ix = nx), deuterium/hydrogen only plasma Z=1
    # check!! if rad_grid from ftxOut.txt has the same convention with SOLPS - python
    Za = 1
    bb, nx, ny = b2fextract('bb',b2fgmtry_file)
    qz, _, _ = b2fextract('qz',b2fgmtry_file)
    pitch = bb[:,:,0]/bb[:,:,3] # bx/bb same as pbs[:,:,0]/qz[:,:,2]/sx - Jeremy
    cosine_Bin = pitch*qz[:,:,1] # pitch * cos t (2D incident angle (on x, y plane)) gives 3D cosine of Bin
    angle_Bin = 57.2958*np.arccos(cosine_Bin) # in degree (180/pi ~ 57.2958 converts unit: rad -> deg)

    te, _, _ = b2fextract('te',b2fstate_file)
    ti, _, _ = b2fextract('ti',b2fstate_file)
    te = te/ec # convert to eV
    ti = ti/ec # convert to eV

    E0_mean = []
    Ain_mean = []
    #pitch_test = []
    #qz_test = []
    angle_Bin_rad_grid = []
    print(rad_grid)
    print(type(rad_grid))
    print(type(rad_grid[0]))
    #rad_grid = range(1,ny) # to test along whole target
    for i in range(0, len(rad_grid)):
        #pitch_test.append(pitch[nx,int(rad_grid[i])])
        #qz_test.append(qz[nx,int(rad_grid[i]),1])
        angle_Bin_rad_grid.append(angle_Bin[nx,int(rad_grid[i])])
        E0_mean.append(3*Za*te[nx,int(rad_grid[i])] + 2*ti[nx,int(rad_grid[i])]) # sum of sheath acceleration and ion kinetic energy
    print(E0_mean)
    print(angle_Bin_rad_grid)
    #print(pitch_test)
    #print(qz_test)





# Sampling of trajectories to get Ein,Ain
# check eirmod_reflec.f


# calculate R_F_TRIM of given sample with Ein, Ain
# also 2D interpolation and extrapolation based on eirmod_reflec.f
# still working on this (20230213)
def get_R_F_TRIM(Ein,Ain,TRIM_file):
    lines = read_lines(TRIM_file)
    for i in range(0, len(lines)):
        # Find closest incident angle among TRIM angle table
        TRIM_angle = np.array([0,30,45,60,70,80,85])
        idx = (np.abs(TRIM_angle-Ain)).argmin()
        angle = TRIM_angle[idx]
    #interpolate or find closest value? how SOLPS deals with it?
        # below number Z, mass? differ by [projectile_on_target] files e.g., 1.  2.00 for D_on_W
        if (len(lines[i].split())>2 and lines[i].split()[0] == '1.' and lines[i].split()[1] == '2.00' and lines[i].split()[5] == angle):
            print (lines[i])
            temp=float(lines[i].split()[6])
            Ein.append(temp_Ein)
            RN.append(temp)

    RN = np.array(RN)
    #np.savetxt(filename+'_'+angle+'_RN_E_mean.txt',np.column_stack([Ein, RN, E_mean]), fmt = ['%.3E','%.3E','%.3E'])
    return


################################## WORKING> ###################################



#Find "*** 6a. General data for reflection model" in input.dat, and update it with FTX output
def input_dat_update(inputdat_file, RECYCF, RECYCT, Twall, rad_grid_name):
    K_to_eV = 8.6173e-05 # boltzmann constant / electron charge to convert K to eV
    # string list to numpy array
    Twall = np.asarray(Twall, dtype = np.float64)
    RECYCF = np.asarray(RECYCF, dtype = np.float64)
    RECYCT = np.asarray(RECYCT, dtype = np.float64)
    Twall = -Twall*K_to_eV # - sign for angular distribution convention in input.dat

    # format numbers in '%.5E' for input.dat e.g., 1.00000E+00
    Twall = ["%.5E" % elem  for elem in Twall] # negative so already 12 digits
    RECYCF = [" %.5E" % elem for elem in RECYCF] # non-negative so need prefix ' ' to complete 12 digits
    RECYCT = [" %.5E" % elem for elem in RECYCT] # non-negative so need prefix ' ' to complete 12 digits
    print('Twall = ', Twall)
    print('RECYCF = ', RECYCF)
    print('RECYCT = ', RECYCT)

    # replace relevant lines in input.dat block 6a
    # find block *** 6a. General data for reflection model
    lines = read_lines(inputdat_file)
    line_block_6a = find_line(lines,'General data for reflection model')
    line_block_7 = find_line(lines, 'Data for primary sources')
    #print(lines[line_block_6a])
    #print(lines[line_block_7])

    surf_index = 0
    for surfmod_name in rad_grid_name: 
        print(surfmod_name)
        # find SURFMOD_name block within block 6a (before block 7)
        line_SURFMOD_def = find_line(lines[line_block_6a:line_block_7], surfmod_name)

        #print before replacement
        #print(lines[line_block_6a+line_SURFMOD_def]) # SURFMOD_name line
        #print(lines[line_block_6a+line_SURFMOD_def+1]) # 1st line ILREF ILSPT ISRS ISRC
        print(lines[line_block_6a+line_SURFMOD_def+2]) # 2nd line ZNML EWALL EWBIN TRANSP(1,N) ...
        print(lines[line_block_6a+line_SURFMOD_def+3]) # 3rd line RECYCF RECYCT RECPRM EXPPL ...

        # change Twall, RECYCF, RECYCT - replace lines in 12 digits input.dat format
        # EWALL: 2nd line, 12:23, RECYCF: 3rd line 0:11, RECYCT: 3rd line 12:23
        lines[line_block_6a+line_SURFMOD_def+2] = replace_str_index(lines[line_block_6a+line_SURFMOD_def+2],12,23,Twall[surf_index]) # replace EWALL (2nd line)
        lines[line_block_6a+line_SURFMOD_def+3] = replace_str_index(lines[line_block_6a+line_SURFMOD_def+3],0,23,RECYCF[surf_index]+RECYCT[surf_index]) # replace RECYCF, RECYCT (3rd line)

        # increase surfmod counter
        surf_index += 1

        # print after replacement of Twall, RECYCF, RECYCT
        print(lines[line_block_6a+line_SURFMOD_def+2])
        print(lines[line_block_6a+line_SURFMOD_def+3])

# write changes to new input.dat (test_inputdat.txt)
    with open('test_inputdat.txt', 'w') as file:
        file.writelines(lines)
    file.close()




### RUN SCRIPTS ###
#AL rad_grid, Twall, R_FT, R_X, R_FTX = get_FTX_out(FTX_file)
#AL RECYCT = R_FTX
#RECYCF = R_FT # just temperarily set as R_FT. will be replaced to the scaling factor later
#AL RECYCF = calc_RECYCF(inputdat_file, fort44_file, rad_grid_name, R_FT)

#AL print('RECYCF',RECYCF)
# To be developed
#RECYCF = scale_factor_RECYCF(R_FT, rad_grid, b2fstate_file) # scale factor R_FT/R_F_TRIM

#AL input_dat_update(inputdat_file, RECYCF, RECYCT, Twall, rad_grid_name)
#AL get_Te_Ti_Bin(rad_grid, b2fstate_file, b2fgmtry_file)

