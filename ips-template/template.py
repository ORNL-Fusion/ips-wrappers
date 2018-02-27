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
        zip_ref.extract(self.services.get_config_param('INPUT_1'))
        zip_ref.extract(self.services.get_config_param('INPUT_2'))
        zip_ref.close()

#  Replace parameters in the extracted state files.
        namelist_example_file = OMFITnamelist(self.services.get_config_param('INPUT_1'))

#  Modify with constant.
        namelist_example_file['example_namelist']['array'][1] = 2.0
#  Modify with keyword arguments.
        for key, val in keywords.iteritems():
#  Some times we only want to modify a specific part of an item. Check if the
#  keyword is a dictionary of indicies and values of the index.
            if isinstance(val, dict):
                for k, v in val.iteritems():
                    namelist_example_file['example_namelist'][key][int(k)] = v
#  Other times modify the entire item.
            namelist_example_file['example_namelist'][key] = val
        
        namelist_example_file.save()

#  Some codes need symbolic links.
        os.symlink(self.services.get_config_param('INPUT_2'), 'link_name')

#-------------------------------------------------------------------------------
#
#  template Component step method. This runs code the component is wrapped
#  around.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        print('template: step')

#  Launch a task. Task launches are non blocking so save the task id. Arguments
#  can be passed in from the config file or passed in the call to the method.
        self.running_tasks['task1'] = self.services.launch_task(self.NPROC,
                                                                self.services.get_working_dir(),
                                                                self.TEMPLATE_EXE,
                                                                *self.TEMPLATE_ARGS.split(' '),
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
