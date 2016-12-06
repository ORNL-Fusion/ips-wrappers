      module b2mod_version
      character*10 :: newversion='01.001.029'
      character*10, save ::
     1 b2fstati_version = '00.000.000',
     2 b2fstate_version = '00.000.000',
     2 b2frates_version = '00.000.000',
     3 b2fmovie_version = '00.000.000',
     3 b2ftrack_version = '00.000.000',
     4 plasmastate_version = '00.000.000',
     4 b2fplasma_version = '00.000.000',
     5 b2fgmtry_version = '00.000.000',
     6 b2fparam_version = '00.000.000',
     7 boundary_version = '00.000.000',
     8 neutrals_version = '00.000.000'
      
      end module b2mod_version 
cxpb Version 1.1 corresponds to adding a version number and the UADIA array
cxpb to the 'b2fstat*' files. Implemented 20.01.2000.
cxpb Version 1.1 also adds the 'periodic_bc' variable to the geometry file.
cxpb Finalised 14.07.2000.
cdpc Version 01.001.001 corresponds to adding b2stbr_sna,smo,she,shi,sch to
cdpc b2mod_sources
cdpc Version 01.001.002 corresponds to adding the anomalous transport 
cdpc (D,chi...) quantities to b2fplasma [2001.07.18]
cdpc Version 01.001.003 corresponds to adding the heat sources  
cdpc from b2sihs to b2fplasma [2002.08.07]
cxpb Version 01.001.004 corresponds to adding the parallel momentum sources
cxpb from atomic processes and the total energy flows to b2fplasma [2002.08.13]
cdpc Version 01.001.005: added saving of fch_32 and fch_52 to b2fstat[ie]
cdpc Version 01.001.006: added saving of b2npmo_sma[fgv] to b2fplasma
cdpc Version 01.001.007: added saving of csigin & cvsa_cl to b2fplasma
cdpc Version 01.001.008: added saving of fmo to b2fplasma
cxpb Version 01.001.009: added saving of kinrgy to b2fstat[ie]
cxpb Version 01.001.010: added saving of wadia and vaecrb to b2fplasma
cdpc Version 01.001.011: added saving of b2stbc_s?? to b2fplasma
cdpc Version 01.001.012: added saving of facdrift and fac_ExB to b2fplasma
cdpc Version 01.001.013: added saving of ext sources
cdpc Version 01.001.014: added saving of additional current arrays
cdpc Version 01.001.015: added saving of fac_vis
cxpb Version 01.001.016: added rtlrd line radiation rate array to b2frates
cdpc Version 01.001.017: added neutral flux limiter diagnostic arrays
cxpb Version 01.001.018: added simulation time to b2fstate (per ANK request)
cxpb Version 01.001.019: added rtlbr bremsstrahlung rate array to b2frates
cxpb Version 01.001.020: added viscosity flux limiter diagnostic array
cxpb Version 01.001.021: added electron and ion heat flux limit arrays
cxpb Version 01.001.022: added qc array to b2fgmtry and b2fplasma
cxp                      added fch_p parallel current array to b2fstat[ie]
cxpb Version 01.001.023: added sclrtio, sclstra arrays to b2fplasma
cxpb Version 01.001.024: added rtlza, rtlz2, rtlpt and rtlpi arrays to b2frates
cxpb                     added saving of ne to b2fstat[ie]
cxpb                     new species data format: zamin, zamax, zn, am
cxpb                     added version number to b2fparam and b2fpardf
cxpb                     added nstra to b2fparam
cxpb                     made ppout a permanent fixture of b2frates
cxpb Version 01.001.025: added pbs computed from basis mesh contributions to b2fgmtry
cxpb Version 01.001.026: added sne, b2stbr_sne, b2stbc_sne and ext_sne to b2fplasma
cxpb Version 01.001.027: added saving of b2stbm_s?? to b2fplasma
cxpb Version 01.001.028: new style for boundary and neutrals namelist to allow for
cxpb                     full vacuum vessel grid coverage
cxpb Version 01.001.029: split the thermo-electric coefficient into a classical and
cxpb                     and an anomalous piece
cxpb                 
