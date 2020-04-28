#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  Example IPS wrapper for FTridyn init component. This wapper only generates
#     input for FTridyn and copies input files to the plasma state
#
#-------------------------------------------------------------------------------

from component import Component
import shutil
import os
import ftMPI
#-------------------------------------------------------------------------------
#
#  FTridyn init Component Constructor
#
#-------------------------------------------------------------------------------
class ftridynInit(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

#-------------------------------------------------------------------------------
#
#  FTridyn init Component init method. This method creates dummy files in the
#    plasma state that can then be updated in the step method.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        print('ftridyn_init: init nothing done this one')
        #Get input file names from config file
        current_ftridyn_namelist = self.services.get_config_param('FTRIDYN_INPUT_FILE')
        #split filenames into a list of strings
        file_list = current_ftridyn_namelist.split()
        #loop over file names and create dummy files in ftridynInit work area
        for index in range(len(file_list)):
            open(file_list[index], 'a').close()

        #Get output file names from config file
        outfile_param = self.services.get_config_param('FTRIDYN_OUTPUT_FILE')
        #split filenames into a list of strings
        outfile_list = outfile_param.split()
        #loop over file names and create dummy files in ftridynInit work area
        for index in range(len(outfile_list)):
            open(outfile_list[index], 'a').close()

        #update plasma state from relevant files in ftridynInit work area
        self.services.update_plasma_state()
        #self.services.stage_input_files(self.INPUT_FILES)

#-------------------------------------------------------------------------------
#
#  FTridyn init Component step method. This runs generateInput.py .
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        print('ftridyn_init: step')

        #call generateInput.py
        dict = {'target':self.TARGET,
                'beam':self.BEAM,
                'nHistories':int(self.N_HISTORIES),
                'incident_energy':float(self.INCIDENT_ENERGY),
                'incident_angle':float(self.ANGLE),
                'data_file':self.DATA_FILE}
        name1 = ''
        name2 = ''
        if(len(dict['beam'])>1):
            name1 = dict['beam']
        else:
            name1 = dict['beam']+'_'
        if(len(dict['target'])>1):
            name2 = dict['target']
        else:
            name2 = '_'+dict['target']

        ftMPI.beam_and_target(name1+name2,dict['beam'],dict['target'],sim_number=1,number_histories=dict['nHistories'], incident_energy=dict['incident_energy'],depth=200.0,incident_angle=dict['incident_angle'],fluence=1.0E-16,dataFile=dict['data_file']) 
        shutil.copyfile(name1+name2+'0001.IN', self.COPY_FILES)
        #get name of FTridyn input file from config file to copy newly generated files to
        current_ftridyn_namelist = self.services.get_config_param('FTRIDYN_INPUT_FILE')
        #this may be more than one file, not sure yet - need to learn more about FTridyn I/O
        from_file_list = self.COPY_FILES.split()
        file_list = current_ftridyn_namelist.split()

        #copy newly generated files to names specified in config file
        for index in range(0,1): #range(len(file_list)): this may need to be changed
            print(('copying ', from_file_list[index], ' to ', file_list[index]))
        #    shutil.copyfile(from_file_list[index], file_list[index])

        #update plasma state files with relevant files from ftridynInit work directory
        self.services.update_plasma_state()
#-------------------------------------------------------------------------------
#
#  FTridyn init Component finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('ftridyn_init: finalize')
