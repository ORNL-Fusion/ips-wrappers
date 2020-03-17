
import numpy as np
import sys
#import matplotlib.pyplot as plt
 
def xolotlToLay(tridynFile='last_TRIDYN.dat', totalDepth=0.0, logFile=None):

  if logFile  is not None:
    print(('redirect tridynPlotting output of to:', logFile))
    outF = open(logFile, "a")
    sys.stdout = outF
    
  else:
    print ('no log file defined in tridynPlotting')
    print ('print output to default sys.stdout')
    
  print(' ')
  print('xolotlToLay:')

  DNS_W_A3_0 = 0.06306
 
  data = np.loadtxt(tridynFile)
  depthnm = data[:,0]
  dx = depthnm[2] - depthnm[1]
  
  DNS_He_nm3 = data[:,1]
  DNS_D_nm3 = data[:,2]
  DNS_T_nm3 = data[:,3]

  DNS_He_A3 = DNS_He_nm3 / 1000.0
  DNS_D_A3 = DNS_D_nm3 / 1000.0
  DNS_T_A3 = DNS_T_nm3 / 1000.0

  C_Total = DNS_W_A3_0 + DNS_He_A3 + DNS_D_A3 + DNS_T_A3
  C_W_A3 = DNS_W_A3_0 / C_Total
  C_He_A3 = DNS_He_A3 / C_Total
  C_D_A3 = DNS_D_A3 / C_Total
  C_T_A3 = DNS_T_A3 / C_Total
  C_prj_A3 = 0.0*C_W_A3  #projectile
 
  #plt.plot(depthnm*10.0,C_He_A3)
  #plt.plot(depthnm*10.0,C_D_A3)
  #plt.plot(depthnm*10.0,C_T_A3)
  #plt.plot(depthnm*10.0,C_W_A3)


#if totalDepth>0.0, add pure W layers until 'totalDepth'  
  if (totalDepth>0.0):
    min_number_lines = int(totalDepth / dx)
    if np.size(C_W_A3)<min_number_lines:
      C_W_A3 = np.hstack((C_W_A3, np.ones(min_number_lines - np.size(C_W_A3))))
      C_He_A3 = np.hstack((C_He_A3, np.zeros(min_number_lines - np.size(C_He_A3))))
      C_D_A3 = np.hstack((C_D_A3, np.zeros(min_number_lines - np.size(C_D_A3))))
      C_T_A3 = np.hstack((C_T_A3, np.zeros(min_number_lines - np.size(C_T_A3))))
      C_prj_A3 = np.hstack((C_prj_A3, np.zeros(min_number_lines - np.size(C_prj_A3))))
    #end if 

  output = open('He_W0001.LAY','w+')
 
  for i in range(np.size(C_W_A3)):
    if(np.isnan(C_W_A3[i]) or np.isnan(C_He_A3[i]) or np.isnan(C_D_A3[i]) or np.isnan(C_T_A3[i])):
      C_W_A3[i] = 0.0
      C_He_A3[i] = 0.0
      C_D_A3[i] = 0.0
      C_T_A3[i] = 0.0
      C_prj_A3[i] = 0.0
    print(C_prj_A3[i],C_W_A3[i],C_He_A3[i],C_D_A3[i],C_T_A3[i],file=output)
  #end for

  nlines=np.size(C_W_A3)
  
  #test in transition to python3:
  print('from translate_xolotl_to_ftridyn, nlines', nlines)
  sys.stdout.flush()
  
  return nlines

#end def main

if __name__ == '__main__':
  xolotlToLay()
