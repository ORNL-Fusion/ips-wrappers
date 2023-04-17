#!/usr/bin/env python

import numpy as np
import sys

def b2fextract(variable,filename, option=None): # e.g. variable = ua, filename='b2fstate'

#    import numpy as np
#    import sys

    is_body=False
    datastr=[]

# Load dimensions

    b2f=open(filename,'r')
    lines = b2f.readlines() # read lines
    for i in range(0, len(lines)):
        if ((lines[i].split()[0] == '*cf:' and lines[i].split()[3] == 'nx,ny,ns') or (lines[i].split()[0] == '*cf:' and lines[i].split()[3] == 'nx,ny')):
            nx=int(lines[i+1].split()[0])
            ny=int(lines[i+1].split()[1])
            #ns=int(lines[i+1].split()[2])
            #print ('nx',nx,'ny',ny,'ns',ns)
            break

# Find a block contains variable

    b2f=open(filename,'r')
    for line in b2f: #.readline(): # line scan from the b2fstate
        #print (line.split()[0])
        #print (line.split()[1])
        #print (line.split()[2])
        #print (line.split()[3])
        if (line.split()[0] == '*cf:' and line.split()[3] == variable):
        # and (int(line.split()[2]) > 0)):
            nentry=int(line.split()[2]) # dimension of the variable
            is_body = True
            ns = int(nentry/(nx+2)/(ny+2))
            print (line)

            #if (variable == 'te' or variable == 'ne' ):
            #    ns=1;

            continue # data splitting and save start from the next line
        if is_body:
            temp=line.split()
            datastr=datastr+temp #  split data by () and append

            if len(datastr) == nentry: break
    
    data = np.array(datastr,dtype=float)
    #data=map(float,datastr) # convert datastr into integer array
    #np.savetxt(variable+'.txt',data)
    if option != 'no_reshape':
        data=np.array(data).reshape(nx+2,ny+2,ns,order='F') # convert data into numpy array and reshape

    #np.save(variable+".npy", data)
    #sys.stdout.write('Output saved as\' '+variable+'.npy\' file')
#    print ('Output saved as\'',variable,'.npy\' file')
    return data

# add executable main part:
# e.g. $python b2fextract.py te b2fstate

#if __name__ == '__main__':
#    variable = sys.argv[1]
#    filename = sys.argv[2]
#    b2fextract(variable,filename)
