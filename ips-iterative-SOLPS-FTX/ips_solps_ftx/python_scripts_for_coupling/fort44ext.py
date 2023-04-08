#!/usr/bin/env python
import math
import numpy as np
import sys

def fort44ext(variable,filename, option=None): # e.g. variable = ua, filename='b2fstate'

#    import numpy as np
#    import sys

    is_body=False
    datastr=[]

# Load dimensions

    b2f=open(filename,'r')
    lines = b2f.readlines() # read lines
    nx = int(lines[0].split()[0])
    ny = int(lines[0].split()[1])
#    print('nx',nx,'ny',ny)
    nentry = 0 # intitialize
#    nx=pol_eirene
#    ny=rad_eirene
    ns=1
#    nentry = 3240# Find a block contains variable

    b2f=open(filename,'r')
    for line in b2f: #.readline(): # line scan from the b2fstate
        if (line.split()[0] == '*eirene' and line.split()[3] == variable): 
            is_body = True
            nentry = int(line.split()[6])
            print (nentry)
            print (line)

            continue # data splitting and save start from the next line
        if is_body:
            temp=line.split()
            datastr=datastr+temp #  split data by () and append

            if len(datastr) == nentry: break

    data = np.array(datastr,dtype=float) # convert datastr into integer array
#    np.savetxt(variable+'.txt',data)
    #data=np.array(data).reshape(nx,ny,ns,order='F') # convert data into numpy array and reshape

    if option != 'no_reshape':
        if nentry > 0 and nentry % (nx * ny) == 0:
            data = np.array(data).reshape(nx, ny, int(nentry / (nx * ny)), order='F')
        elif nentry == 0 or math.isnan(nentry):
            print('Error: variable not found in fort.44 file')
        else:
            pass




    #np.save(variable+".npy", data)
    #sys.stdout.write('Output saved as\' '+variable+'.npy\' file')
#    print 'Output saved as\'',variable,'.npy\' file'
    return data

# add executable main part:
# e.g. $python b2fextract.py te b2fstate
#if __name__ == '__main__':
#    variable = sys.argv[1]
#    filename = sys.argv[2]
#    pol_eirene = sys.argv[3]
#    rad_eirene = sys.argv[4]
#    fort44ext(variable,filename,pol_eirene,rad_eirene)
#    fort44ext(variable,filename)

