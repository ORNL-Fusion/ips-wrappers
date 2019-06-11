#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for TEMPLATE component. This template shows an example of how to
#  set up an ips component.
#
#-------------------------------------------------------------------------------

from component import Component
import os
from omfit.classes.omfit_namelist import OMFITnamelist
import zipfile

#-------------------------------------------------------------------------------
#
#  TEMPLATE Class
#
#-------------------------------------------------------------------------------
class template(Component):
    
#-------------------------------------------------------------------------------
#
#  TEMPLATE Init Component Constructor
#
#-------------------------------------------------------------------------------
    def __init__(self, services, config):
        print('template: Construct')
        Component.__init__(self, services, config)
        self.running_tasks = {}

#-------------------------------------------------------------------------------
#
#  template Component init method. This method prepairs the input files.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0, **keywords):
        print('template: init')

        #  Stage the plasma state files.
        self.services.stage_plasma_state()

        #  Unzip files from the plasma state.
        zip_ref = zipfile.ZipFile(self.services.get_config_param('CURRENT_TEMPLATE_STATE'), 'r')
        zip_ref.extract(self.services.get_config_param('TEMPLATE_NAMELIST_INPUT'))
        zip_ref.extract(self.services.get_config_param('TEMPLATE_DAT_INPUT'))
        zip_ref.close()

        #  Modify entries of a namelist input file using OMFIT.
        namelist_example_file = OMFITnamelist(self.services.get_config_param('TEMPLATE_NAMELIST_INPUT'))

        #  Modify using data specified within this wrapper.
        namelist_example_file['example_namelist']['array'][1] = 2.0

        #  Modify using data passed as keyword arguments.
        for key, val in keywords.iteritems():
            
            #  Some times we only want to modify a specific part of an item. Check if the
            #  keyword is a dictionary of indicies and values of the index.
            
            if isinstance(val, dict):
                for k, v in val.iteritems():
                    namelist_example_file['example_namelist'][key][int(k)] = v
            
            #  Other times modify the entire item.
            namelist_example_file['example_namelist'][key] = val
        
        # Modify using data from ips.config file
        
            # to be implemented
        
        # Modify using data from a PPPL plasma state file
        
            # to be implemented
            
        # Modify using data from an OMAS state file
        
            # to be implemented
        
        # Write the resulting namelist file
        namelist_example_file.save()

        #  Some codes need symbolic links. This will throw an expection if the link
        #  allready exists.
        if os.path.isfile('link_name'):
            os.remove('link_name')
        os.symlink(self.services.get_config_param('TEMPLATE_DAT_INPUT'), 'link_name')
            
#-------------------------------------------------------------------------------
#
#  template Component step method. This runs code the component is wrapped
#  around.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0, **keywords):
        print('template: step')

        #  Arguments can be passed in from the config file or passed in the call to the
        #  step method. All definitions from the config file come in as a string.
        #  Keywords maynot not exist so check if the key exists first.
        template_args = [ self.TEMPLATE_ARGS ]
        if 'TEMPLATE_ARGS' in keywords:
            template_args.extend(keywords['TEMPLATE_ARGS'])

        #  Launch a task. Task launches are non blocking so save the task id. Depending
        #  on the system, install locations and binary names can be different
        self.running_tasks['task1'] = self.services.launch_task(self.NPROC,
                                                                self.services.get_working_dir(),
                                                                self.TEMPLATE_EXE,
                                                                *template_args,
                                                                logfile = 'template.log')
                                                                
        #  Wait for the task to complete.
        self.services.wait_task(self.running_tasks['task1'])
        del self.running_tasks['task1']

#-------------------------------------------------------------------------------
#
#  template Component finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('template: finalize')
        for value in self.running_tasks:
            self.services.wait_task(value)

