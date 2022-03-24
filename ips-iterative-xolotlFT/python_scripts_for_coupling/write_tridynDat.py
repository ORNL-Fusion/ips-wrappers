#!/usr/bin/env python
# write tridyn.dat file
#format determined by version of Xolotl


import os
import sys

def write_tridynDat(outFile='tridyn.dat', tridynDat_model=1, plasmaSpecies=['He','W','D','T'], timeFolder='t0.0', maxRangeXolotl=[0.0, 0.0, 0.0, 0.0], fluxFraction=[0.0, 0.0, 0.0, 0.0], rYield=[1.0, 1.0, 1.0, 1.0], xp_parameters={}  ): #INPUTS HERE
    print(' ')
    print('from write_tridynDat, called with input')
    print('\t outFile =', outFile)
    print('\t tridynDat_model =',tridynDat_model )
    print('\t plasmaSpecies =',plasmaSpecies )
    print('\t timeFolder =', timeFolder)
    print('\t maxRangeXolotl =', maxRangeXolotl)
    print('\t fluxFraction =', fluxFraction)
    print('\t rYield =', rYield)
    print('\t xp_parameters =', xp_parameters)
    sys.stdout.flush()
    
    ## New attempt at checking which species should be written into tridyn.dat:
    ##  for tridynDat_model==2,
    ##    - first check for flux fraction ;
    ##             if 0, skip writing into tridyn.dat because there's no information on impact energy and angle (no FT)
    ##                   regardless of network: lines missing for species that exist in the network does NOT cause any issues (just nothing implanted)
    ##             if >0, check network (#netParam = nHe nD nT maxVSize nInterstitials) -->
    ##    - use the network param values when available,
    ##    - assume the standard (species included if fluxFraction>0, not included if fluxFraction<0) if the network isn't given explicitely
    ##      I.e., now we merge two checks: it works in cases where more species passed by GITR, but not handled by Xolotl; or if fluxFraction=0.
    ##  for tridynDat_model==1, include all species:
    ##   - zero if fluxFraction==0, or if netParam == 0 (when available)

    if os.path.exists(outFile):
        os.remove(outFile)
    combinedFile = open(outFile, 'a')
    
    #tridynDat_model==2: tridyn.dat format for (e.g.) pulsed, UQ & tempGrid executables of Xolotl:
    #                    header: prj cluster_size (1-refl)
    #                    only for all prj present in plasma & networkFile
    if (tridynDat_model==2):
        for i in range(len(plasmaSpecies)):
            prj=plasmaSpecies[i]
            ft_output_profile_temp_prj=timeFolder+'/'+outFile+'_'+prj
            profile=open(ft_output_profile_temp_prj, "r")
            tridynString=profile.read().rstrip('\n')
            combinedTridynString=str(tridynString)+str(maxRangeXolotl[i])
            profile.close()
            print(('for {0}, fraction in plasma = {1} , and reflection = {2} '.format(prj,fluxFraction[i], rYield[i])))
            print(('\t effective fraction (in plasma * (1-reflection)) = {} '.format(fluxFraction[i]*(1-rYield[i]))))
            sys.stdout.flush()
            
            if (fluxFraction[i] > 0):
                print('\t Write tridyn.dat line for ', prj, ', in new tridyn.dat format (model ', str(tridynDat_model),')')
                # if not W, then the name in the tridyn.dat line is the same as prj
                # if He,  check He's position in netParam, i.e., index i=0
                if prj=='He':
                    if ('netParam' in xp_parameters):
                        if (xp_parameters['netParam'][i]==0):
                            print('\t Xolotl netowrk exists for ' , prj, 'given in plasmaSpecies, but entry in netParam is zero ; will skip in tridyn.dat')
                        else:
                            print('\t For ' , prj , 'netparam = ' , xp_parameters['netParam'][i] , ' is used in Xolotl ; write line for ', prj , ' in tridyn.dat')
                            combinedFile.write("%s %s %s\n" %(prj,str(1),str(fluxFraction[i]*(1-rYield[i]))))
                            combinedFile.write("%s\n" %(combinedTridynString))
                    else:
                        print('\t WARNING:',prj,' exist in plasma but netparam not given in Xolotl ; write line for in tridyn.dat')
                        print('\t \t this might give an ERROR if species isnt part of Xolotls network')
                        combinedFile.write("%s %s %s\n" %(prj,str(1),str(fluxFraction[i]*(1-rYield[i]))))
                        combinedFile.write("%s\n" %(combinedTridynString))
                ##if D or T,  check position in netParam, i.e., index i=2,3 -> i-1 = 1 or 2 (no W in netParam)
                elif prj=='D' or prj=='T':
                    if ('netParam' in xp_parameters):
                        if (xp_parameters['netParam'][i-1]==0):
                            print('\t Xolotl network exists for ' , prj, 'given in plasmaSpecies, but entry in netParam is zero ; will skip in tridyn.dat')
                        else:
                            print('\t For ' , prj , 'netparam = ' , xp_parameters['netParam'][i-1] , ' is used in Xolotl ; write line for ', prj , ' in tridyn.dat')
                            combinedFile.write("%s %s %s\n" %(prj,str(1),str(fluxFraction[i]*(1-rYield[i]))))
                            combinedFile.write("%s\n" %(combinedTridynString))
                    else:
                        print('\t WARNING:',prj,' exist in plasma but netparam not given in Xolotl ; write line for in tridyn.dat')
                        print('\t \t this might give an ERROR if species isnt part of Xolotls network')
                        combinedFile.write("%s %s %s\n" %(prj,str(1),str(fluxFraction[i]*(1-rYield[i]))))
                        combinedFile.write("%s\n" %(combinedTridynString))
                elif prj=='W':
                    if ('netParam' in xp_parameters):
                        if (xp_parameters['netParam'][4]==0):
                            print('\t Xolotl netowrk exists for ' , prj, 'given in plasmaSpecies, but entry in netParam is zero ; will skip in tridyn.dat')
                        else:
                            print('\t For ' , prj , 'netparam = ' , xp_parameters['netParam'][4] , ' is used in Xolotl ; write line for ', prj , ' in tridyn.dat')
                            combinedFile.write("%s %s %s\n" %('I',str(1),str(fluxFraction[i]*(1-rYield[i]))))
                            combinedFile.write("%s\n" %(combinedTridynString))
                    else:
                        print('\t WARNING:',prj,' exist in plasma but netparam not given in Xolotl ; write line for in tridyn.dat')
                        print('\t \t this might give an ERROR if species isnt part of Xolotls network')
                        combinedFile.write("%s %s %s\n" %('I',str(1),str(fluxFraction[i]*(1-rYield[i]))))
                        combinedFile.write("%s\n" %(combinedTridynString))
                    
                else:
                    print('\t WARNING: species ', prj, 'cannot be handled by Xolotl yet.')
                    print('\t \t it has been used so far (for spY, etc), but will skip writing into tridyn.dat')

            elif (fluxFraction[i] == 0):
                print('\t Using the new tridyn.dat format (model = ', str(tridynDat_model), '), flux fraction for ', prj, 'is zero')
                print('\t \t for now, skip writing anything, even if prj exists in network (no checks in place)')
                print('\t \t Xolotl will run, with no ', prj, ' implanted')
                sys.stdout.flush()
        sys.stdout.flush()
                
    #tridynDat_model==1: tridyn.dat format for (e.g.) master executable of Xolotl:
    #                    header: (1-refl)
    #                    for all prj in He W D T
    elif (tridynDat_model==1):
        xolSpecies=['He','W','D','T']
        for i in range(len(xolSpecies)):
            prj=xolSpecies[i]
            print('\t Write tridyn.dat line for ', prj, ', in original tridyn.dat format (model ', str(tridynDat_model),')')
            ft_output_profile_temp_prj=timeFolder+'/'+outFile+'_'+prj
            if os.path.exists(ft_output_profile_temp_prj):
                profile=open(ft_output_profile_temp_prj, "r")
                tridynString=profile.read().rstrip('\n')
                combinedTridynString=str(tridynString)+str(maxRangeXolotl[i])
                print(('for {0}, fraction in plasma = {1} , and reflection = {2} '.format(prj,fluxFraction[i], rYield[i])))
                print(('\t effective fraction (in plasma * (1-reflection)) = {} '.format(fluxFraction[i]*(1-rYield[i]))))
                profile.close()
            else:
                print('no ', ft_output_profile_temp_prj, ' found; set tridyn.dat values to zero')
                combinedTridynString='0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n'
                sys.stdout.flush()

            # if He,  check He's position in netParam, i.e., index i=0
            if (prj=='He'):
                if ('netParam' in xp_parameters):
                    if (xp_parameters['netParam'][i]==0):
                        print('\t For ' , prj , 'netparam = ' , xp_parameters['netParam'][i] , ' --> set all entries to 0.0 in tridyn.dat')
                        combinedFile.write("0.0\n")
                        combinedFile.write("0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n")
                    else:
                        print('\t For ' , prj , 'netparam = ' , xp_parameters['netParam'][i] , ' in Xolotl ; write line in tridyn.dat')
                        combinedFile.write("%s\n" %(str(fluxFraction[i]*(1-rYield[i]))))
                        combinedFile.write("%s\n" %(combinedTridynString))
                else:
                    print('\t WARNING: netparam not given in Xolotl ; write line for ', prj , ' in tridyn.dat')
                    print('\t \t this might give an ERROR if species isnt part of Xolotls network')
                    combinedFile.write("%s\n" %(str(fluxFraction[i]*(1-rYield[i]))))
                    combinedFile.write("%s\n" %(combinedTridynString))
            # if W,  check W's position in netParam, i.e., index i=4
            elif (prj=='W'):
                if ('netParam' in xp_parameters):
                    if (xp_parameters['netParam'][4]==0):
                        print('\t For ' , prj , 'netparam = ' , xp_parameters['netParam'][4] , ' --> set all entries to 0.0 in tridyn.dat')
                        combinedFile.write("0.0\n")
                        combinedFile.write("0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n")
                    else:
                        print('\t For ' , prj , 'netparam = ' , xp_parameters['netParam'][4] , ' in Xolotl ; write line in tridyn.dat')
                        combinedFile.write("%s\n" %(str(fluxFraction[i]*(1-rYield[i]))))
                        combinedFile.write("%s\n" %(combinedTridynString))
                else:
                    print('\t WARNING: netparam not given in Xolotl ; write line for ', prj , ' in tridyn.dat')
                    print('\t \t this might give an ERROR if species isnt part of Xolotls network')
                    combinedFile.write("%s\n" %(str(fluxFraction[i]*(1-rYield[i]))))
                    combinedFile.write("%s\n" %(combinedTridynString))
            ##if D or T,  check position in netParam, i.e., index i=2,3 -> i-1 = 1 or 2
            elif (prj=='D') or (prj=='T'):
                if ('netParam' in xp_parameters):
                    if (xp_parameters['netParam'][i-1]==0):
                        print('\t For ' , prj , 'netparam = ' , xp_parameters['netParam'][i-1] , ' --> set all entries to 0.0 in tridyn.dat')
                        combinedFile.write("0.0\n")
                        combinedFile.write("0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n")
                    else:
                        print('\t For ' , prj , 'netparam = ' , xp_parameters['netParam'][i-1] , ' in Xolotl ; write line in tridyn.dat')
                        combinedFile.write("%s\n" %(str(fluxFraction[i]*(1-rYield[i]))))
                        combinedFile.write("%s\n" %(combinedTridynString))
                else:
                    print('\t WARNING: netparam not given in Xolotl ; write line for ', prj , ' in tridyn.dat')
                    print('\t \t this might give an ERROR if species isnt part of Xolotls network')
                    combinedFile.write("%s\n" %(str(fluxFraction[i]*(1-rYield[i]))))
                    combinedFile.write("%s\n" %(combinedTridynString))
            else:
                print('\t WARNING: species ', prj, ' cannot be handled by Xolotl yet.')
                print('\t \t it has been used so far (for spY, etc), but will skip writing into tridyn.dat')
            sys.stdout.flush()
    else:
        print('\t WARNING: tridynDat_model ', str(tridynDat_model), ' not recognized.')
        print('\t print in original (master) format, with all zeros')
        combinedFile.write("0.0\n")
        combinedFile.write("0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n")
    combinedFile.close()

    print('... done with write_tridynDat')

    sys.stdout.flush()
