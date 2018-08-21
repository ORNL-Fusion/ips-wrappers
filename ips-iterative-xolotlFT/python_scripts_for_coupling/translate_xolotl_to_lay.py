from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
 
def xolotlToLay(tridynFile='TRIDYN_last.dat'):
  DNS_W_A3_0 = 0.06306
 
  data = np.loadtxt(tridynFile)
  depthnm = data[:,0]
  DNS_He_nm3 = data[:,1]
  dx = depthnm[1] - depthnm[2]
  DNS_He_A3 = DNS_He_nm3 / 1000.0
 
  C_Total = DNS_W_A3_0 + DNS_He_A3
  C_W_A3 = DNS_W_A3_0 / C_Total
  C_He_A3 = DNS_He_A3 / C_Total
 
  #plt.plot(depthnm*10.0,C_He_A3)
  #plt.plot(depthnm*10.0,C_W_A3)
 
  output = open('He_W0001.LAY','w+')
 
  for i in range(np.size(C_W_A3)):
    if(np.isnan(C_W_A3[i]) or np.isnan(C_He_A3[i])):
      C_W_A3[i] = 0.0
      C_He_A3[i] = 0.0
    print(C_He_A3[i],C_W_A3[i],file=output)
  #end for
#end def main

  nlines=len(depthnm)

  return nlines

if __name__ == '__main__':
  xolotlToLay()
