#! /usr/bin/env python

import os
import f90nml

from  component import Component

class cheasepy_driver(Component):
      def __init__(self, services, config):
          Component.__init__(self, services, config)
          print('Created %s' % (self.__class__))

      def init(self, timeStamp = 0.0):
          
          print('LAUNCH: cheasepy_driver_init')

          self.inputfpath = getattr(self,'INPUT_DIR',  '')
          self.inputfname = getattr(self,'INPUT_FILES','inchease0') 
          if os.path.isfile('%s/%s' %(self.inputfpath,self.inputfname)):
             inchease = f90nml.read('%s/%s' %(self.inputfpath,self.inputfname))
          else:
             inchease = {}

          if 'sources'  in inchease:
             srcVals     = dict(inchease['sources'])
          else:
             srcVals     = {}
          if 'namelist' in inchease:
             nlVals = dict(inchease['namelist'])
             namelistVals = {}
             for ikey in nlVals:
                 namelistVals[ikey.upper()] = nlVals[ikey]
          else:
             namelistVal = {}

          if 'rhomesh_src' not in srcVals:
             inptval = getattr(self,'rhomesh_src',      'eqdsk')
             srcVals['rhomesh_src']   = int(float(inptval)) if inptval.isnumeric() else inptval.strip()
          if 'current_src' not in srcVals:
             inptval = getattr(self,'current_src',      'eqdsk')
             srcVals['current_src']   = int(float(inptval)) if inptval.isnumeric() else inptval.strip()
          if 'boundary_src' not in srcVals:
             inptval = getattr(self,'boundary_src',      'asis')
             srcVals['boundary_src']  = int(float(inptval)) if inptval.isnumeric() else inptval.strip()
          if 'pressure_src' not in srcVals:
             inptval = getattr(self,'pressure_src',     'eqdsk')
             srcVals['pressure_src']  = int(float(inptval)) if inptval.isnumeric() else inptval.strip()
          if 'eprofiles_src' not in srcVals:
             inptval = getattr(self,'eprofiles_src', 'profiles')
             srcVals['eprofiles_src'] = int(float(inptval)) if inptval.isnumeric() else inptval.strip()
          if 'iprofiles_src' not in srcVals:
             inptval = getattr(self,'iprofiles_src', 'profiles')
             srcVals['iprofiles_src'] = int(float(inptval)) if inptval.isnumeric() else inptval.strip()
          if 'inputpath' not in srcVals:
             inptval = getattr(self,'inputpath',         '')
             srcVals['inputpath']     = int(float(inptval)) if inptval.isnumeric() else inptval.strip()
          if 'gfname' not in srcVals:
             inptval = getattr(self,'gfname',            '')
             srcVals['gfname']        = int(float(inptval)) if inptval.isnumeric() else inptval.strip()
          if 'pfname' not in srcVals:
             inptval = getattr(self,'pfname',            '')
             srcVals['pfname']        = int(float(inptval)) if inptval.isnumeric() else inptval.strip()
          if 'expeqfname' not in srcVals:
             inptval = getattr(self,'expeqfname',        '')
             srcVals['expeqfname']    = int(float(inptval)) if inptval.isnumeric() else inptval.strip()
          if 'cheasefname' not in srcVals:
             inptval = getattr(self,'cheasefname',       '')
             srcVals['cheasefname']   = int(float(inptval)) if inptval.isnumeric() else inptval.strip()
          if 'exptnzfname' not in srcVals:
             inptval = getattr(self,'exptnzfname',       '')
             srcVals['exptnzfname']   = int(float(inptval)) if inptval.isnumeric() else inptval.strip()
          if 'iterdbfname' not in srcVals:
             inptval = getattr(self,'iterdbfname',       '')
             srcVals['iterdbfname']   = int(float(inptval)) if inptval.isnumeric() else inptval.strip()

          if 'NS' not in namelistVals:
             namelistVals['NS']        = int(float(getattr(self,'NS',       256)))
          if 'NT' not in namelistVals:
             namelistVals['NT']        = int(float(getattr(self,'NT',       256)))
          if 'NISO' not in namelistVals:
             namelistVals['NISO']      = int(float(getattr(self,'NISO',     256)))
          if 'NPSI' not in namelistVals:
             namelistVals['NPSI']      = int(float(getattr(self,'NPSI',    1024)))
          if 'NCHI' not in namelistVals:
             namelistVals['NCHI']      = int(float(getattr(self,'NCHI',    1024)))
          if 'NRBOX' not in namelistVals:
             namelistVals['NRBOX']     = int(float(getattr(self,'NRBOX',     60)))
          if 'NZBOX' not in namelistVals:
             namelistVals['NZBOX']     = int(float(getattr(self,'NZBOX',     60)))
          if 'RELAX' not in namelistVals:
             namelistVals['RELAX']     =     float(getattr(self,'RELAX',      0)) 
          if 'NSTTP' not in namelistVals:
             namelistVals['NSTTP']     = int(float(getattr(self,'NSTTP',      1)))
          if 'NPROPT' not in namelistVals:
             namelistVals['NPROPT']    = int(float(getattr(self,'NPROPT',     3)))
          if 'NPPFUN' not in namelistVals:
             namelistVals['NPPFUN']    = int(float(getattr(self,'NPPFUN',     8)))
          if 'NEQDSK' not in namelistVals:
             namelistVals['NEQDSK']    = int(float(getattr(self,'NEQDSK',     1)))
          if 'NFUNRHO' not in namelistVals:
             namelistVals['NFUNRHO']   = int(float(getattr(self,'NFUNRHO',    0)))
          if 'TENSBND' not in namelistVals:
             namelistVals['TENSBND']   = int(float(getattr(self,'TENSBND',    0)))
          if 'TNESPROF' not in namelistVals:
             namelistVals['TENSPROF']  = int(float(getattr(self,'TENSPROF',   0)))
          if 'NRHOMESH' not in namelistVals:
             namelistVals['NRHOMESH']  = int(float(getattr(self,'NRHOMESH',   0)))
          if 'cocos_in' not in namelistVals:
             namelistVals['cocos_in']  = int(float(getattr(self,'cocos_in',   2)))
          if 'cocos_out' not in namelistVals:
             namelistVals['cocos_out'] = int(float(getattr(self,'cocos_out', 12)))

          self.services.stage_input_files(self.INPUT_FILES)
          try:
              worker_comp = self.services.get_port('WORKER')
          except Exception:
              self.services.exception('Error accessing worker component')
              raise
          self.services.call(worker_comp, 'init', 0.0, sources=srcVals, namelist=namelistVals)
          print('FINISH: cheasepy_driver_int')

          return

      def step(self, timeStamp=0.0):

          print('LAUNCH: cheasepy_driver_step')
          try:
              worker_comp = self.services.get_port('WORKER')
          except Exception:
              self.services.exception('Error accessing worker component')
              raise
          self.services.call(worker_comp, 'step',     0.0)
          self.services.call(worker_comp, 'finalize', 0.0)
          print('FINISH: cheasepy_driver_step')

          return

      def finalize(self, timeStamp=0.0):
          return

