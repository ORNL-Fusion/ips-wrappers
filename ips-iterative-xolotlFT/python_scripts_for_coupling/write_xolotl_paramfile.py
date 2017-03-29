#!/usr/bin/env python
#=======================================================================================
# write_xoltol_paramfile.py
# writes the parameter file for xolotl (outputfile), replacing variable values (tmax, network file, sputtering yield) in the input file (infile)
#=======================================================================================

import os
import subprocess


##  v1
############### ORIGINAL BASH SCRIPT ###########

#infile=`echo "param.template"`;
#outfile=`echo "param.txt"`;
#dt=0.3
#netfile=`echo "networkxolotl"`
#sputtyield=0.1

#loadtxt=('infile')

#eval "sed -e 's/-start_stop tmax/-start_stop $dt/' -e 's/-ts_final_time tmax/-ts_final_time $dt/' -e 's/networkFile=xolotlnetworkfile/networkFile=$netfile/' -e 's/sputtering=sputteringyield/sputtering=$sputtyield/ '<param.template >param.txt "

########### END ORIGINAL BASH SCRIPT (v1) ###########



##  v2
########### ORIGINAL PYTHON SCRIPT, MODIFYING A TEMPLATE FILE ################

def writeXolotlParameterFile_fromTemplate(infile="params.txt", outfile="params.txt",

                                          start_stop=0.01,
                                          
                                          ts_final_time=0.01,
                                          ts_max_snes_failures=-1,
                                          ts_max_steps=1000000,
                                          ts_exact_final_time="matchstep",
                                          ts_adapt_dt_max=1.0e-5,
                                          ts_monitor=True,
                                          
                                          fieldsplit_0_pc_type="sor",
                                          fieldsplit_1_pc_type="redundant",

                                          pc_fieldsplit_detect_coupling=True,
                                          pc_type="fieldsplit",
                                          
                                          vizHandler="dummy",
                                          flux=4.0e4,
                                          networkFile="networkInit.h5",
                                          material="TRIDYN",
                                          dimensions=1,
                                          perfHandler="dummy",
                                          startTemp=773,
                                          process="reaction advec modifiedTM attenuation diff movingSurface",
                                          sputtering=0.000129
                                          ):
   tmp="temp.txt"
#   ftmp=open("temp.txt", "w")                                                                                                                                                                
   if (infile==outfile):
      os.rename(infile, tmp)
      infile=tmp

# prepare petscline:                                                                                                                                                                          

   #change of value in parameters                                                                                                                                                             
   petscArgString=" -e 's/-start_stop [^ ]*/-start_stop %f/'   -e 's/-ts_final_time [^ ]*/-ts_final_time %f/'   -e 's/-ts_max_snes_failures [^ ]*/-ts_max_snes_failures %d/'   -e 's/-ts_max_steps [^ ]*/-ts_max_steps %d/'   -e 's/-ts_exact_final_time [^ ]*/-ts_exact_final_time %s/'   -e 's/-ts_adapt_dt_max [^ ]*/-ts_adapt_dt_max %e/'   -e 's/-fieldsplit_0_pc_type [^ ]*/-fieldsplit_0_pc_type %s/'   -e 's/-fieldsplit_1_pc_type [^ ]*/-fieldsplit_1_pc_type %s/'   -e 's/-pc_type [^ ]*/-pc_type %s/' "  % (start_stop, ts_final_time, ts_max_snes_failures, ts_max_steps, ts_exact_final_time, ts_adapt_dt_max, fieldsplit_0_pc_type, fieldsplit_1_pc_type, pc_type)

   #include (or not) parameters without values that exist in template file

#   if (ts_monitor==True):
#      petscArgString=petscArgString+"   -e 's/-ts_monitor/-ts_monitor/'"
#   else:
#      petscArgString=petscArgString+"   -e 's/-ts_monitor/ /'"

#   if (pc_fieldsplit_detect_coupling==True):
#      petscArgString=petscArgString+"   -e 's/-pc_fieldsplit_detect_coupling/ -pc_fieldsplit_detect_coupling/'"
#   else:
#      petscArgString=petscArgString+"   -e 's/-pc_fieldsplit_detect_coupling/ /'"

   #prepare sed line                                                                                                                                                                          
   petscArgSedString="sed "+ petscArgString + "< %s > %s" %(infile , outfile)

   #run sed line for Petsc                                                                                                                                                                    

   subprocess.call([petscArgSedString], shell=True)


   #other input parameters

   os.rename(outfile, tmp)
   paramSedString="sed    -e 's/vizHandler=[^ ]*/vizHandler=%s/'    -e 's/flux=[^ ]*/flux=%e/'    -e 's/networkFile=.*$/networkFile=%s/'    -e 's/material=[^ ]*/material=%s/'    -e 's/dimensions=[^ ]*/dimensions=%d/'    -e 's/perfHandler=[^ ]*/perfHandler=%s/'    -e 's/startTemp=[^ ]*/startTemp=%f/'   -e 's/sputtering=[^ ]*/sputtering=%f/'    -e 's/process=.*$/process=%s/' < %s > %s"   % (vizHandler, flux, networkFile, material, dimensions, perfHandler, startTemp, sputtering, process, tmp, outfile)

   #print " sedline call parameters: %s " %(paramSedString)                                                                                                                                   
   subprocess.call([paramSedString], shell=True)

   os.remove(tmp)

   return






#   sedstring="sed -e 's/-start_stop tmax/-start_stop %f/' -e 's/-ts_final_time tmax/-ts_final_time %f/' -e 's/networkFile=xolotlnetworkfile/networkFile=%s/' -e 's/sputtering=\sputteringyield/sputtering=%f/ '< %s > %s " % (tmax, tmax, networkFile, sputtering, infile , outfile)

#   print "\n running the following sed command: %s " %(sedstring)                                                                                                                            
#   subprocess.call([sedstring], shell=True)


##################### END ORIGINAL PYTHON SCRIPT (v2)  ####################



##  v3
########## NEW PYTHON SCRIPT: MODIFY PARAMTER FILE CREATED BY XOLOTL ##############

def writeXolotlParameterFile_fromPreprocessor(infile="params.txt", outfile="params.txt",

                         tridyn=True,
                         helium_retention=True,

                         start_stop=0.01,

                         ts_final_time=0.01,
                         ts_max_snes_failures=-1,
                         ts_max_steps=1000000,
                         ts_exact_final_time="matchstep",
                         ts_adapt_dt_max=1.0e-5,
                         ts_monitor=True,

                         fieldsplit_0_pc_type="sor",
                         fieldsplit_1_pc_type="redundant",

                         pc_fieldsplit_detect_coupling=True,
                         pc_type="fieldsplit",

                         vizHandler="dummy",
                         flux=4.0e4,
                         networkFile="networkInit.h5",
                         material="TRIDYN",
                         dimensions=1,
                         perfHandler="dummy",
                         startTemp=773,
                         process="reaction advec modifiedTM attenuation diff movingSurface",
                         regularGridInline=True,
                         regularGrid="no",
                         voidPortionInline=True,
                         voidPortion=40.0,
                         sputteringInline=True,
                         sputtering=0.000129
                         ):


   "write the parameter file for xolotl replacing tmax, networkfile and sputtering yield values in the input file"
   "for tests, run: python write_xolotl_paramfile.py"

   tmp="temp.txt"
#   ftmp=open("temp.txt", "w")

   if (infile==outfile):
      os.rename(infile, tmp)
      infile=tmp


# prepare petscline:

   #change of value in parameters
   petscArgString=" -e 's/-ts_dt 1.0e-12/-start_stop %f/'   -e 's/-ts_final_time 1.0/-ts_final_time %f/'   -e 's/-ts_max_snes_failures 200/-ts_max_snes_failures %d/'   -e 's/-ts_max_steps 100/-ts_max_steps %d/'   -e 's/-ts_exact_final_time stepover/-ts_exact_final_time %s/'   -e 's/-ts_adapt_dt_max 1.0e-6/-ts_adapt_dt_max %e/'   -e 's/-fieldsplit_0_pc_type sor/-fieldsplit_0_pc_type %s/'   -e 's/-fieldsplit_1_pc_type redundant/-fieldsplit_1_pc_type %s/'   -e 's/-pc_type fieldsplit/-pc_type %s/' "  % (start_stop, ts_final_time, ts_max_snes_failures, ts_max_steps, ts_exact_final_time, ts_adapt_dt_max, fieldsplit_0_pc_type, fieldsplit_1_pc_type, pc_type)

   #include (or not) parameters without values that exist in template file
   if (ts_monitor==True):
      petscArgString=petscArgString+"   -e 's/-ts_monitor/-ts_monitor/'"
   else:
      petscArgString=petscArgString+"   -e 's/-ts_monitor/ /'"

   if (pc_fieldsplit_detect_coupling==True):
      petscArgString=petscArgString+"   -e 's/-pc_fieldsplit_detect_coupling/ -pc_fieldsplit_detect_coupling/'"
   else:
      petscArgString=petscArgString+"   -e 's/-pc_fieldsplit_detect_coupling/ /'"

   #prepare sed line 
   petscArgSedString="sed "+ petscArgString + "< %s > %s" %(infile , outfile)

   #run sed line for Petsc
   
   subprocess.call([petscArgSedString], shell=True)

   #append parameters not present in template file to Petsc line 

   if (tridyn==True):
      os.rename(outfile, tmp)
      petscArgAppendSed1="sed '/petscArgs=/s/$/  -tridyn/' < %s > %s " %(tmp, outfile) 
      subprocess.call([petscArgAppendSed1], shell=True)

   if (helium_retention==True):
      os.rename(outfile, tmp)
      petscArgAppendSed2="sed '/^petscArgs=/ s/$/  -helium_retention/'  < %s > %s " %(tmp, outfile)
      subprocess.call([petscArgAppendSed2], shell=True)

   #other input parameters

   os.rename(outfile, tmp)
   paramSedString="sed    -e 's/vizHandler=dummy/vizHandler=%s/'    -e 's/flux=4.0e7/flux=%e/'    -e 's/networkFile=networkInit.h5/networkFile=%s/'    -e 's/material=W100/material=%s/'    -e 's/dimensions=1/dimensions=%d/'    -e 's/perfHandler=std/perfHandler=%s/'    -e 's/startTemp=1000/startTemp=%f/'    -e 's/process=reaction diff advec/process=%s/' < %s > %s"   % (vizHandler, flux, networkFile, material, dimensions, perfHandler, startTemp, process, tmp, outfile)

   #print " sedline call parameters: %s " %(paramSedString)
   subprocess.call([paramSedString], shell=True)

   #append other input parameters that do not exist in preprocessors param file

   f = open(outfile, "a")

   if (regularGridInline==True):
      regularGridLine="regularGrid=%s\n" %(regularGrid)
      f.write(regularGridLine)

   if (voidPortionInline==True):
      voidPortionLine="voidPortion=%f\n" %(voidPortion) 
      f.write(voidPortionLine)

   if (sputteringInline==True):
      sputteringLine="sputtering=%f\n" %(sputtering)
      f.write(sputteringLine)

   f.close

   os.remove(tmp)

   return 

################# END OF NEW PYTHON SCRIPT (v3) ####################

if __name__ == '__main__':

#   infile="params.txt" 
#   outfile="params.txt"
#   networkFile="networkInit.testing.txt"
#   sputtering=0.129
   
#   writexolotlparamfile(infile=infile, outfile=outfile, networkFile=networkFile, sputtering=sputtering)

   import shutil

   writeXolotlParameterFile_fromPreprocessor()

   shutil.copyfile("params.txt", "params1.txt")

   writeXolotlParameterFile_fromTemplate(start_stop=0.2,ts_final_time=0.2,networkFile="xolotlStop.h5",sputtering=0.1)

   shutil.copyfile("params.txt", "params2.txt")
