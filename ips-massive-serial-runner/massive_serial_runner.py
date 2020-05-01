#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for the Massive Serial Runner component.
#
#-------------------------------------------------------------------------------

from component import Component
from utilities import ScreenWriter
from utilities import ZipState
import os
import json
import numpy
import netCDF4
from omfit.classes.omfit_eqdsk import OMFITeqdsk

#-------------------------------------------------------------------------------
#
#  Linear integration.
#
#-------------------------------------------------------------------------------
def lin_integral(x, y):
    return numpy.add.reduce(0.0625*(x[1:] + x[:-1])*(y[1:] + y[:-1])*(x[i:] - x[:-1])

#-------------------------------------------------------------------------------
#
#  Dictionary of constant values.
#
#-------------------------------------------------------------------------------
constants = {
    'i_core'   : 50,
    'nus_edge' : 0.0,
    'mu0'      : math.pi*4.0E-7
}

#-------------------------------------------------------------------------------
#
#  Massive Serial Runner Component Constructor.
#
#  The parameter is a dictionary of fastran outputs and their prefix codes.
#  Prefix refers to the method for populating that value. Some prefix codes are
#  actualy post processing routines and contain dependancies. Name is the
#  internal name used in the code. This can differ from the key. For post
#  processed parameters, the depends lists the dependancy of the parameter needs
#  to be loaded before this value can be computed.
#
#-------------------------------------------------------------------------------
class massive_serial_runner(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)

        self.parameter = {
            'peak'        : {'prefix' : self.aeqdisk,        'name' : 'peak',        'depends' : []                                                     },
            'betap'       : {'prefix' : self.aeqdisk,        'name' : 'betap',       'depends' : []                                                     },
            'betan_efit'  : {'prefix' : self.aeqdisk,        'name' : 'betan',       'depends' : []                                                     },
            'q95'         : {'prefix' : 'unknown',           'name' : 'q95',         'depends' : []                                                     },
            'li'          : {'prefix' : 'unknown',           'name' : 'li',          'depends' : []                                                     },
            'qmin'        : {'prefix' : 'unknown',           'name' : 'qmin',        'depends' : []                                                     },
            'rho_qmin'    : {'prefix' : 'unknown',           'name' : 'rho_qmin',    'depends' : []                                                     },
            'ne'          : {'prefix' : self.fastran,        'name' : 'ne',          'depends' : []                                                     },
            'ni'          : {'prefix' : self.fastran,        'name' : 'ni',          'depends' : []                                                     },
            'te'          : {'prefix' : self.fastran,        'name' : 'te',          'depends' : []                                                     },
            'ti'          : {'prefix' : self.fastran,        'name' : 'ti',          'depends' : []                                                     },
            'q'           : {'prefix' : self.fastran,        'name' : 'q',           'depends' : []                                                     },
            'rmajor'      : {'prefix' : self.fastran,        'name' : 'rmajor',      'depends' : []                                                     },
            'aminor'      : {'prefix' : self.fastran,        'name' : 'aminor',      'depends' : []                                                     },
            'zeff'        : {'prefix' : self.fastran,        'name' : 'zeff',        'depends' : []                                                     },
            'kappa'       : {'prefix' : self.fastran,        'name' : 'delta',       'depends' : []                                                     },
            'delta'       : {'prefix' : self.fastran,        'name' : 'kappa',       'depends' : []                                                     },
            'ip'          : {'prefix' : self.fastran,        'name' : 'ip',          'depends' : []                                                     },
            'bt'          : {'prefix' : self.fastran,        'name' : 'b0',          'depends' : []                                                     },
            'r'           : {'prefix' : self.fastran,        'name' : 'r0',          'depends' : []                                                     },
            'a'           : {'prefix' : self.fastran,        'name' : 'a0',          'depends' : []                                                     },
            'we'          : {'prefix' : self.fastran,        'name' : 'we',          'depends' : []                                                     },
            'wi'          : {'prefix' : self.fastran,        'name' : 'wi',          'depends' : []                                                     },
            'wb'          : {'prefix' : self.fastran,        'name' : 'wb',          'depends' : []                                                     },
            'betan'       : {'prefix' : self.fastran,        'name' : 'betan',       'depends' : []                                                     },
            'nebar'       : {'prefix' : self.fastran,        'name' : 'nebar',       'depends' : []                                                     },
            'tea'         : {'prefix' : self.fastran,        'name' : 'tea',         'depends' : []                                                     },
            'tia'         : {'prefix' : self.fastran,        'name' : 'tia',         'depends' : []                                                     },
            'taue'        : {'prefix' : self.fastran,        'name' : 'taue',        'depends' : []                                                     },
            'taui'        : {'prefix' : self.fastran,        'name' : 'taui',        'depends' : []                                                     },
            'tauth'       : {'prefix' : self.fastran,        'name' : 'tauth',       'depends' : []                                                     },
            'tautot'      : {'prefix' : self.fastran,        'name' : 'tautot',      'depends' : []                                                     },
            'tautot'      : {'prefix' : self.fastran,        'name' : 'tautot',      'depends' : []                                                     },
            'tau89'       : {'prefix' : self.fastran,        'name' : 'tau89',       'depends' : []                                                     },
            'tau98'       : {'prefix' : self.fastran,        'name' : 'tau98',       'depends' : []                                                     },
            'pnbe'        : {'prefix' : self.fastran,        'name' : 'pnbe',        'depends' : []                                                     },
            'pnbi'        : {'prefix' : self.fastran,        'name' : 'pnbi',        'depends' : []                                                     },
            'prfe'        : {'prefix' : self.fastran,        'name' : 'prfe',        'depends' : []                                                     },
            'prfi'        : {'prefix' : self.fastran,        'name' : 'prfi',        'depends' : []                                                     },
            'prad'        : {'prefix' : self.fastran,        'name' : 'prad',        'depends' : []                                                     },
            'pei'         : {'prefix' : self.fastran,        'name' : 'pei',         'depends' : []                                                     },
            'poh'         : {'prefix' : self.fastran,        'name' : 'poh',         'depends' : []                                                     },
            'ibs'         : {'prefix' : self.fastran,        'name' : 'ibs',         'depends' : []                                                     },
            'inb'         : {'prefix' : self.fastran,        'name' : 'inb',         'depends' : []                                                     },
            'irf'         : {'prefix' : self.fastran,        'name' : 'irf',         'depends' : []                                                     },
            'pfuse_equiv' : {'prefix' : self.fastran,        'name' : 'pfuse_equiv', 'depends' : []                                                     },
            'pfusi_equiv' : {'prefix' : self.fastran,        'name' : 'pfusi_equiv', 'depends' : []                                                     },
            'rho'         : {'prefix' : self.fastran,        'name' : 'rho',         'depends' : []                                                     },
            'rhob'        : {'prefix' : self.fastran,        'name' : 'rhob',        'depends' : []                                                     },
            'volp'        : {'prefix' : self.fastran,        'name' : 'volp',        'depends' : []                                                     },
            'pe_rf'       : {'prefix' : self.fastran,        'name' : 'pe_rf',       'depends' : []                                                     },
            'area'        : {'prefix' : self.fastran,        'name' : 'area',        'depends' : []                                                     },
            'pfus0'       : {'prefix' : self.outone,         'name' : 'pfus0',       'depends' : []                                                     },
            'tauR'        : {'prefix' : self.get_tauR,       'name' : 'tauR',        'depends' : ['sigma', 'rho', 'rhob']                               },
            'sigma'       : {'prefix' : self.sigma_hirshman, 'name' : 'sigma',       'depends' : ['ne', 'te', 'aminor', 'rmajor', 'zeff', 'q']          },
            'neavg'       : {'prefix' : self.get_neavg,      'name' : 'sigma',       'depends' : ['volp', 'rhob', 'ne']                                 },
            'rho_ec'      : {'prefix' : self.get_rho_ec,     'name' : 'rho_ec',      'depends' : ['rho', 'pe_rf']                                       },
            'q0'          : {'prefix' : self.get_q0,         'name' : 'q0',          'depends' : ['q']                                                  },
            'nus_core'    : {'prefix' : self.nues,           'name' : 'nus_core',    'depends' : ['rmajor', 'aminor', 'ne', 'te', 'zeff', 'q', 'i_core']},
            'kappa1'      : {'prefix' : self.get_kappa1,     'name' : 'kappa1',      'depends' : ['kappa']                                              },
            'delta1'      : {'prefix' : self.get_delta1,     'name' : 'delta1',      'depends' : ['delta']                                              },
            'area1'       : {'prefix' : self.get_area1,      'name' : 'area1',       'depends' : ['r', 'a', 'kappa', 'delta']                           },
            'pinj'        : {'prefix' : self.get_pinj,       'name' : 'pinj',        'depends' : ['prfe', 'prfi', 'pnbe', 'pnbi']                       },
            'pinj_ei'     : {'prefix' : self.get_pinj_ei,    'name' : 'pinj_ei',     'depends' : ['prfe', 'pnbe', 'pinj']                               },
            'pfus'        : {'prefix' : self.get_pfus,       'name' : 'pfus',        'depends' : ['pfuse_equiv', 'pfusi_equiv']                         },
            'qdt'         : {'prefix' : self.get_qdt,        'name' : 'qdt',         'depends' : ['pfus', 'pinj']                                       },
            'wmhd'        : {'prefix' : self.get_wmhd,       'name' : 'wmhd',        'depends' : ['we', 'wi', 'wb']                                     },
            'fbs'         : {'prefix' : self.get_fbs,        'name' : 'fbs',         'depends' : ['ibs', 'ip']                                          },
            'fnb'         : {'prefix' : self.get_fnb,        'name' : 'fnb',         'depends' : ['inb', 'ip']                                          },
            'frf'         : {'prefix' : self.get_frf,        'name' : 'frf',         'depends' : ['irf', 'ip']                                          },
            'fni'         : {'prefix' : self.get_fni,        'name' : 'fni',         'depends' : ['fbs', 'fnb', 'frf']                                  },
            'h98'         : {'prefix' : self.get_h98,        'name' : 'h98',         'depends' : ['tauth', 'tau98']                                     },
            'fgw'         : {'prefix' : self.get_fgw,        'name' : 'fgw',         'depends' : ['nebar', 'ip', 'a']                                   },
            'fgw_ped'     : {'prefix' : self.get_fgw_ped,    'name' : 'fgw_ped',     'depends' : ['neped', 'ip', 'a']                                   },
            'pmhd'        : {'prefix' : self.get_pmhd,       'name' : 'pmhd',        'depends' : ['b', 'ip', 'a', 'betan', 'mu0']                       },
            'nus_edge'    : {'prefix' : self.constant,       'name' : 'nus_edge',    'depends' : []                                                     },
            'i_core'      : {'prefix' : self.constant,       'name' : 'i_core',      'depends' : []                                                     },
            'mu0'         : {'prefix' : self.constant,       'name' : 'mu0',         'depends' : []                                                     }
        }

        self.open_fastran = []
        self.open_aeqdisk = []

#-------------------------------------------------------------------------------
#
#  Massive Serial Runner Component init method. This method prepairs the state.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'massive_serial_runner: init')

#  Get config filenames.
        if (timeStamp == 0.0):
            self.current_state = self.services.get_config_param('CURRENT_MSR_STATE')
            self.current_batch = self.services.get_config_param('CURRENT_MSR_BATCH')
            self.model_config = self.services.get_config_param('MODEL_CONFIG')
            self.fastran_config = self.services.get_config_param('FASTRAN_CONFIG')
            self.shot_number = '00000' # FIXME: this is a placeholder.

#  Stage state.
        self.services.stage_state()

#  Unzip files from the state. Use mode a so files an be read and written to.
        self.zip_ref = ZipState.ZipState(self.current_state, 'a')
        self.zip_ref.extract('inscan')

#-------------------------------------------------------------------------------
#
#  Massive Serial Runner Component step method. Not used.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'massive_serial_runner: step')

        flags = self.zip_ref.get_state()

#  Run the massive serial workflow.
        if 'state' in flags and flags['state'] == 'needs_update':
            with open(current_batch, 'r') as json_ref:
                keys = json.load(json_ref).keys()

            task_wait = self.services.launch_task(MAX(self.NPROC, len(keys)),
                                                  self.services.get_working_dir(),
                                                  self.services.MASSIVE_SERIAL_EXE,
                                                  'inscan',
                                                  self.fastran_config,
                                                  MAX(self.NPROC, len(keys)),
                                                  1)

#  Collect results of the workflow.
            if (self.services.wait_task(task_wait)):
                self.services.error('massive_serial_runner: step failed.')

#  Append output data to batch. Need to know what data to read.

#-------------------------------------------------------------------------------
#
#  Massive Serial Runner Component finalize method. This cleans up
#  afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'massive_serial_runner: finalize')

#-------------------------------------------------------------------------------
#
#  Load data parameter. Data parameters are listed in the parameters dictionary.
#  for a given parameter, start by loading the dependencies. Then based on the
#  prefix value load the data from the appropriate source. To avoid file io
#  overhead, a file will be first checked it is opened. Opened files will be
#  added to a list to cache the file reference. This updates the data base with
#  the newly loaded parameter.
#
#  FIXME: Not sure how I should update for the mutilple directories so start with just one.
#
#-------------------------------------------------------------------------------
    def load_parameter(self, directory, parameter, data_base, batch_size):

#  First check that the parameter is not already loaded.
        if parameter in data_base:
            return data_base

#  Next load any dependencies.
        for key in parameters[parameter][depends]:
            data_base = self.load_parameter(directory, key, parameters)

#  Now load the parameter.
        data_base[parameter] = [self.parameters[parameter]['prefix'](self.parameters[parameter]['name'], data_base, i) for i in range(batch_size)]

#  Close the open datasets.


#-------------------------------------------------------------------------------
#
#  Load fastran parameter.
#
#-------------------------------------------------------------------------------
    def fastran(self, parameter, data_base, index):
        if len(self.open_fastran) <= index:
            file_name = 'out{:05d}/f{}.{:05d}'.format(index, self.shot_number, index)
            self.open_fastran.append(netCDF4.Dataset(file_name,'r', format='NETCDF4'))

        return self.open_fastran[index][parameter][:,:]

#-------------------------------------------------------------------------------
#
#  Load a parameter from a EFIT a file.
#
#-------------------------------------------------------------------------------
    def aeqdisk(self, parameter, data_base, index):
        if len(self.open_fastran) <= index:
            file_name = 'out{:05d}/a{}.{:05d}'.format(index, self.shot_number, index)
            self.open_aeqdisk.append(OMFITeqdsk(file_name))

        return self.open_aeqdisk[index][parameter]

#-------------------------------------------------------------------------------
#
#  Load a parameter from a outone file. This file only contains one parameter so
#  there is no need to cache the open file.
#
#-------------------------------------------------------------------------------
    def outone(self, parameter, data_base, index):
        pattern = re.compile("\s*P DD")

        file_name = 'out{:05d}/a{}.{:05d}'.format(index, self.shot_number, index)
        with open(file_name, 'r') as outone_ref:
            for line in outone_ref.readlines():
                if pattern.search(line):
                    return float(line.split()[2][1:])*222.0/1.0E6

#-------------------------------------------------------------------------------
#
#  Load a parameter from a constant.
#
#-------------------------------------------------------------------------------
    def constant(self, parameter, data_base, index):
        return constants[parameter]

#-------------------------------------------------------------------------------
#
#  Load sigma parameters. This is from Hirshman, Nuclear Fusion 21, 1079 (1981)
#
#  NOTE:
#     lnlamda ~ 24.0 - log (sqrt(ne)/te) in N/cm^3, eV
#     nu_star ~ sqrt(2)*R0*q*nu_ee/eps^1.5/vth_e
#
#  These variable names are not discriptive.
#
#-------------------------------------------------------------------------------
    def sigma_hirshman(self, parameter, data_base, index):
        rtor = data_base['rmajor'][index][-1]
        eps = data_base['aminor'][index][-1]/rtor
        ft = 1.0 - (math.sqrt(1.0 - eps*eps)(1.0 - eps)**2)/(1.0 + 1.46*math.sqrt(eps))

        ne = data_base['ne'][index][-1]
        te = data_base['te'][index][-1]
        lnlamda = 15.941 - 0.5*math.log(ne) + math.log(te)

        zeff = data_base['zeff'][index][-1]
        fz1 = (1.0 + 2.966*zeff + 0.753*zeff*zeff)
        fz2 = (1.0 + 1.198*zeff + 0.222*zeff*zeff)

        q = data_base['q'][index][-1]
        sigma_spitzer = 3.0619E2*zeff*fz1/(lnlamda*(te**1.5)*fz2)
        nu_star = 6.9298E-5*lnlamda*ne*rtor*q/(te*te*eps**1.5)

        zeta = 0.58 + 0.2*zeff
        cr = 0.56*(3.0 - zeff)/(zeff*(3.0 + zeff))

        temp = ft/(1.0 + zeta*nu_star)

        return sigma_spitzer*(1.0 - temp)*(1.0 - cr*temp)

#-------------------------------------------------------------------------------
#
#  Load get_tauR parameters.
#
#-------------------------------------------------------------------------------
    def get_tauR(self, parameter, data_base, index):
        rho = data_base['rho'][index][:]
        rhob = data_base['rhob'][index][-1][-1]
        return 0.17*rhob*rhob*lin_integral(rho, sigma)

#-------------------------------------------------------------------------------
#
#  Load get_tauR parameters.
#
#-------------------------------------------------------------------------------
    def get_neavg(self, parameter, data_base, index):
        rhob = data_base['rhob'][index][-1]
        drho = rhob[1] - rhob[0]

        volp = data_base['volp'][index][-1]
        vol = numpy.zeros(len(rhob))
        for i in range(1, len(rhob)):
            vol[i] = vol[i - 1] + drho*(volp[k - 1] + volp[k])*0.5

        ne = data_base['ne'][index]
        return numpy.add.reduce((vol[1:] - vol[:-1])*ne[-1][:-1])/vol[-1]

#-------------------------------------------------------------------------------
#
#  Load get_rho_ec parameters.
#
#-------------------------------------------------------------------------------
    def get_rho_ec(self, parameter, data_base, index):
        rho = data_base['rho'][index][:]
        pec = data_base['pe_rf'][index][-1]
        return rho[nunpy.argmax(pec)]

#-------------------------------------------------------------------------------
#
#  Load get_q0 parameters.
#
#-------------------------------------------------------------------------------
    def get_q0(self, parameter, data_base, index):
        return data_base['q'][index][-1][0]

#-------------------------------------------------------------------------------
#
#  Load nues parameters.
#
#-------------------------------------------------------------------------------
    def nues(self, parameter, data_base, index):
        icore = 50

        R = data_base['rmajor'][index][-1][i_core]
        r = data_base['aminor'][index][-1][i_core]
        ne = data_base['ne'][index][-1][i_core]
        te = data_base['te'][index][-1][i_core]
        zeff = data_base['zeff'][index][-1][i_core]
        q = data_base['q'][index][-1][i_core]

        coulg = 15.9 - 0.5*math.log(ne) + math.log(te)
        nuee = 670.0*coulg*ne/(te**1.5))
        eps = r/R
        vthe = (te**1.5)*1.875E7
        return 1.4*zeff*nuee*R*q*vthe/eps**1.5

#-------------------------------------------------------------------------------
#
#  Load kappa1 parameters.
#
#-------------------------------------------------------------------------------
    def get_kappa1(self, parameter, data_base, index):
        return data_base['kappa'][index][-1][-1]

#-------------------------------------------------------------------------------
#
#  Load delta1 parameters.
#
#-------------------------------------------------------------------------------
    def get_delta1(self, parameter, data_base, index):
        return data_base['delta'][index][-1][-1]

#-------------------------------------------------------------------------------
#
#  Load area1 parameters.
#
#-------------------------------------------------------------------------------
    def get_area1(self, parameter, data_base, index):
        r = data_base['r'][index][-1]
        a = data_base['a'][index][-1]
        kappa = data_base['kappa'][index][-1]
        delta = data_base['delta'][index][-1]
        return 4.0*math.pi*math.pi*a*kappa*(1.0 - 0.151*delta*a/r) - 0.151*2.0*math.pi*math.pi*a*a*kappa*delta/r

#-------------------------------------------------------------------------------
#
#  Load pinj parameters.
#
#-------------------------------------------------------------------------------
    def get_pinj(self, parameter, data_base, index):
        return data_base['prfe'][index] + data_base['prfi'][index] + data_base['pnbe'][index] + data_base['pnbi'][index]

#-------------------------------------------------------------------------------
#
#  Load pinj_ei parameters.
#
#-------------------------------------------------------------------------------
    def get_pinj_ei(self, parameter, data_base, index):
        return data_base['prfe'][index] + data_base['pnbe'][index])/data_base['pinj'][index]

#-------------------------------------------------------------------------------
#
#  Load pfus parameters.
#
#-------------------------------------------------------------------------------
    def get_pfus(self, parameter, data_base, index):
        return 5.0*(data_base['pfuse_equiv'][index]  + data_base['pfusi_equiv'][index] )

#-------------------------------------------------------------------------------
#
#  Load qdt parameters.
#
#-------------------------------------------------------------------------------
    def get_qdt(self, parameter, data_base, index):
        return data_base['pfus'][index] /(data_base['pinj'][index]  - 0.2*data_base['pfus'][index])

#-------------------------------------------------------------------------------
#
#  Load wmhd parameters.
#
#-------------------------------------------------------------------------------
    def get_wmhd(self, parameter, data_base, index):
        return data_base['we'][index] + data_base['wi'][index] + data_base['wb'][index]

#-------------------------------------------------------------------------------
#
#  Load fni parameters.
#
#-------------------------------------------------------------------------------
    def get_fni(self, parameter, data_base, index):
        return data_base['fbs'][index] + data_base['fnb'][index] + data_base['frf'][index]

#-------------------------------------------------------------------------------
#
#  Load fbs parameters.
#
#-------------------------------------------------------------------------------
    def get_fbs(self, parameter, data_base, index):
        return data_base['ibs'][index]/data_base['ip'][index]

#-------------------------------------------------------------------------------
#
#  Load fnb parameters.
#
#-------------------------------------------------------------------------------
    def get_fnb(self, parameter, data_base, index):
        return data_base['inb'][index]/data_base['ip'][index]

#-------------------------------------------------------------------------------
#
#  Load frf parameters.
#
#-------------------------------------------------------------------------------
    def get_frf(self, parameter, data_base, index):
        return data_base['irf'][index]/data_base['ip'][index]

#-------------------------------------------------------------------------------
#
#  Load h98 parameters.
#
#-------------------------------------------------------------------------------
    def get_h98(self, parameter, data_base, index):
        return data_base['tauth'][index]/data_base['tau98'][index]

#-------------------------------------------------------------------------------
#
#  Load fgw parameters.
#
#-------------------------------------------------------------------------------
    def get_fgw(self, parameter, data_base, index):
        return data_base['nebar'][index]*math.pi*data_base['a0'][index]*data_base['a0'][index]/(10.0*data_base['ip'][index])

#-------------------------------------------------------------------------------
#
#  Load fgw_ped parameters.
#
#-------------------------------------------------------------------------------
    def get_fgw_ped(self, parameter, data_base, index):
        return data_base['neped'][index]*math.pi*data_base['a0'][index]*data_base['a0'][index]/(10.0*data_base['ip'][index])

#-------------------------------------------------------------------------------
#
#  Load pmhd parameters.
#
#-------------------------------------------------------------------------------
    def get_pmhd(self, parameter, data_base, index):
        return  data_base['b0'][index]*data_base['ip'][index]*data_base['betan'][index]*1.0e-5/(data_base['a0'][index]*2.0*data_base['mu0'][index])
