#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  Example IPS wrapper for FTridyn init component. This wapper only generates
#     input for FTridyn and copies input files to the plasma state
#
#-------------------------------------------------------------------------------

from ipsframework import Component
import shutil
import os
import sys
import generateInputIPS
import translate_xolotl_to_lay
#-------------------------------------------------------------------------------
#
#  FTridyn init Component Constructor
#
#-------------------------------------------------------------------------------
class ftridynInit(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

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
       
        param_cfg_file = self.services.get_config_param('PARAMETER_CONFIG_FILE')
        open(param_cfg_file,'a').close 
        #Get xolotl file names from config file
        xolotl_param = self.services.get_config_param('CURRENT_XOLOTL')
        #split filenames into a list of strings
        xolotl_list = xolotl_param.split()
        #loop over file names and create dummy files in ftridynInit work area
        for index in range(len(xolotl_list)):
            open(xolotl_list[index], 'a').close()

        #Get output file names from config file
        transfer_files = self.services.get_config_param('XFT_TRANSFER_FILE')
        #split filenames into a list of strings
        transfer_list = transfer_files.split()
        #loop over file names and create dummy files in ftridynInit work area
        for index in range(len(transfer_list)):
            open(transfer_list[index], 'a').close()

        #Get output file names from config file
        other_files = self.services.get_config_param('OTHER_FILES')
        #split filenames into a list of strings
        other_list = other_files.split()
        #loop over file names and create dummy files in ftridynInit work area
        for index in range(len(other_list)):
            open(other_list[index], 'a').close()
        #update plasma state from relevant files in ftridynInit work area
        self.services.update_state()
        #self.services.stage_input_files(self.INPUT_FILES)
        sys.path.append(os.getcwd())
        import parameterConfig
#-------------------------------------------------------------------------------
#
#  FTridyn init Component step method. This runs generateInput.py .
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        print('ftridyn_init: step')
        self.services.stage_state()

        import parameterConfig
        reload(parameterConfig)
        #parameter_config_file = self.services.get_config_param('PARAMETER_CONFIG_FILE')
        print('from parameterConfig file ', parameterConfig.mode)

        if (parameterConfig.mode == 'INIT'):
            print('init mode yes')
            #call generateInput.py
            #os.system(' '.join(['python', self.INPUT_SCRIPT, '-R 1 -s 0']))
            generateInputIPS.main()
        else:
            print('init mode no')
            nDataPts = translate_xolotl_to_lay.xolotlToLay()
            generateInputIPS.main(IQ0=-1,NQX=nDataPts)
        #get name of FTridyn input file from config file to copy newly generated files to
        current_ftridyn_namelist = self.services.get_config_param('FTRIDYN_INPUT_FILE')
        #this may be more than one file, not sure yet - need to learn more about FTridyn I/O
        from_file_list = self.COPY_FILES.split()
        file_list = current_ftridyn_namelist.split()

        #copy newly generated files to names specified in config file
        for index in range(0,1): #range(len(file_list)): this may need to be changed
            print('copying ', from_file_list[index], ' to ', file_list[index])
            shutil.copyfile(from_file_list[index], file_list[index])
        #update plasma state files with relevant files from ftridynInit work directory
        self.services.update_state()
#-------------------------------------------------------------------------------
#
#  FTridyn init Component finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('ftridyn_init: finalize')
