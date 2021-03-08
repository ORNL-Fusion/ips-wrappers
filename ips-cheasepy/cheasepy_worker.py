#! /usr/bin/env python

import os
import f90nml
import cheasepy

from glob      import glob
from component import Component

class cheasepy_worker(Component):
      def __init__(self, services, config):
          Component.__init__(self, services, config)
          print('Created %s' % (self.__class__))

      def init(self, timeStamp=0.0, **kwargs):
          if kwargs:
             for ikwarg in kwargs:
                 if   ikwarg == 'sources':  self.srcVals      = kwargs['sources']
                 elif ikwarg == 'namelist': self.namelistVals = kwargs['namelist']
          self.importedVals = {}
          cheasepy.init_chease_inputs(self.srcVals,self.namelistVals,self.importedVals)
          return

      def step(self, timeStamp=0.0):
          self.services.stage_plasma_state()
          self.state_work_dir  = self.services.get_config_param("STATE_WORK_DIR")
          self.cur_state_file  = self.services.get_config_param("CURRENT_STATE")
          self.cur_eqdsk_file  = self.services.get_config_param("CURRENT_EQDSK")
          self.cur_chease_file = self.services.get_config_param("CURRENT_CHEASE")

          chease_bin  = os.path.join(self.BIN_PATH, self.BIN)
          cwd         = self.services.get_working_dir()

          print("Starting Initial Setup:")
          task_id     = self.services.launch_task(1,cwd,chease_bin,logfile='xchease.log')
          retcode     = self.services.wait_task(task_id)
          if retcode != 0: raise Exception('Error Executing CHEASE')

          os.system('mv EQDSK_COCOS_02_POS.OUT %s' % self.cur_eqdsk_file )
          os.system('mv ogyropsi.dat %s'           % self.cur_chease_file)
          os.system('mv EXPTNZ DIIID_162940_EXPTNZ'                        )
          os.system('mv EXPEQ  DIIID_162940_EXPEQ'                         )
          self.services.update_plasma_state()

          self.inputfpath = getattr(self,'INPUT_DIR',  '')
          self.inputfname = getattr(self,'INPUT_FILES','inchease1') 
          if os.path.isfile('%s/%s' %(self.inputfpath,self.inputfname)):
             inchease = f90nml.read('%s/%s' %(self.inputfpath,self.inputfname))
          else:
             inchease = {}

          if 'sources'  in inchease:
             srcVals     = dict(inchease['sources'])
             for ikey in self.srcVals:
                 if ikey in srcVals: 
                    self.srcVals[ikey] = srcVals[ikey]

          if 'namelist' in inchease:
             namelistVals = dict(inchease['namelist'])
             for ikey in self.namelistVals:
                 if ikey.lower() in namelistVals:
                    self.namelistVals[ikey.upper()] = namelistVals[ikey.lower()]

          self.srcVals['inputpath']     = self.state_work_dir
          self.srcVals['eqdskfname']    = self.cur_eqdsk_file
          self.srcVals['cheasefname']   = self.cur_chease_file
          self.srcVals['profilesfname'] = self.cur_state_file

          self.srcVals['expeqfname']    = "DIIID_162940_EXPEQ"
          self.srcVals['exptnzfname']   = "DIIID_162940_EXPTNZ"

          for i in range(2):
              print("Starting Iteration(%d):" % (i+1))
              cheasepy.init_chease_inputs(self.srcVals,self.namelistVals,self.importedVals)
              task_id     = self.services.launch_task(1,cwd,chease_bin,logfile='xchease.log')
              retcode     = self.services.wait_task(task_id)
              if retcode != 0: raise Exception('Error Executing CHEASE')

              os.system('mv EQDSK_COCOS_02_POS.OUT %s' % self.cur_eqdsk_file)
              os.system('mv ogyropsi.dat %s'           % self.cur_chease_file)
              os.system('mv EXPTNZ DIIID_162940_EXPTNZ'                      )
              os.system('mv EXPEQ  DIIID_162940_EXPEQ'                       )
              self.services.update_plasma_state()
          cheasepy.plot_chease(cheasefpath=self.cur_chease_file,eqdskfpath=self.inputfpath+"/DIIID_162940_EQDSK",
                               expeqfpath='DIIID_162940_EXPEQ',exptnzfpath='DIIID_162940_EXPTNZ',
                               profilesfpath=self.cur_state_file)

          return

#-------------------------------------------------------------------------------
#  CHEASEPY_init init Component finalize method. This cleans up afterwards.
#-------------------------------------------------------------------------------
      def finalize(self, timeStamp=0.0):

          return
