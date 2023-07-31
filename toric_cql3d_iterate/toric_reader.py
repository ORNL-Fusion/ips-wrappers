import numpy as np
import scipy.io as spio
import string
import matplotlib.pyplot as plt

class toric_file():
    def __init__(self,toric_name='fort.9',lean=True):
        #open cql3d netcdf
        #-----------------------------------------------------------------------------
        try:
            toric_nc = spio.netcdf_file(toric_name,'r')
        except:
            print('toric_file initialization failed: could not find ncdf: ',toric_name)
            raise Exception('toric_file initialization failed: could not find ncdf: ',toric_name)

        #read in cdf dimensions
        #-----------------------------------------------------------------------------
        self.n_of_field_comp  = np.copy(toric_nc.dimensions['n_of_field_comp'])
        self.n_of_pol_modes   = np.copy(toric_nc.dimensions['n_of_pol_modes'])
        self.n_of_pol_pts     = np.copy(toric_nc.dimensions['n_of_pol_pts'])
        self.n_of_pol_modvac  = np.copy(toric_nc.dimensions['n_of_pol_modvac'])
        self.n_of_rad_elem    = np.copy(toric_nc.dimensions['n_of_rad_elem'])
        self.n_of_ionspec     = np.copy(toric_nc.dimensions['n_of_ionspec'])
        self.dim2_of_celem    = np.copy(toric_nc.dimensions['dim2_of_celem'])
        self.n_of_rad_pts     = np.copy(toric_nc.dimensions['n_of_rad_pts'])
        self.n_of_mhd_modes   = np.copy(toric_nc.dimensions['n_of_mhd_modes'])
        self.n_of_mhd_rad_pts = np.copy(toric_nc.dimensions['n_of_mhd_rad_pts'])
        self.data_plawall     = np.copy(toric_nc.dimensions['data_plawall'])

        #transcripe some dims to shorthand
        self.ndims = self.n_of_field_comp
        self.nmod  = self.n_of_pol_modes
        self.ntt   = self.n_of_pol_pts
        self.nelm  = self.n_of_rad_pts
        self.nmhd  = self.n_of_mhd_rad_pts

        #read in cdf variables
        #-----------------------------------------------------------------------------
        if(lean==False):
            self.enhcol = np.copy(toric_nc.variables['enhcol'].data)
            self.dnures = np.copy(toric_nc.variables['dnures'].data)
            self.tnures = np.copy(toric_nc.variables['tnures'].data)
            self.ant_length = np.copy(toric_nc.variables['ant_length'].data)
            self.ant_constant = np.copy(toric_nc.variables['ant_constant'].data)
            self.ant_position = np.copy(toric_nc.variables['ant_position'].data)
            self.torus_radius = np.copy(toric_nc.variables['torus_radius'].data)
            self.axis_radius = np.copy(toric_nc.variables['axis_radius'].data)
            self.plasma_radius = np.copy(toric_nc.variables['plasma_radius'].data)
            self.sep_radius = np.copy(toric_nc.variables['sep_radius'].data)
            self.farshield_radius = np.copy(toric_nc.variables['farshield_radius'].data)
            self.antenna_radius = np.copy(toric_nc.variables['antenna_radius'].data)
            self.wall_radius = np.copy(toric_nc.variables['wall_radius'].data)
            self.b_zero = np.copy(toric_nc.variables['b_zero'].data)
            self.b_axis = np.copy(toric_nc.variables['b_axis'].data)
            self.tor_current = np.copy(toric_nc.variables['tor_current'].data)
            self.psi_edge = np.copy(toric_nc.variables['psi_edge'].data)
            self.ppjte = np.copy(toric_nc.variables['ppjte'].data)
            self.ppjti = np.copy(toric_nc.variables['ppjti'].data)
            self.Shfr_shift_axis = np.copy(toric_nc.variables['Shfr_shift_axis'].data)
            self.Shfr_shift_wall = np.copy(toric_nc.variables['Shfr_shift_wall'].data)
            self.ellip_axis = np.copy(toric_nc.variables['ellip_axis'].data)
            self.ellip_edge = np.copy(toric_nc.variables['ellip_edge'].data)
            self.ellip_wall = np.copy(toric_nc.variables['ellip_wall'].data)
            self.triang_edge = np.copy(toric_nc.variables['triang_edge'].data)
            self.triang_wall = np.copy(toric_nc.variables['triang_wall'].data)
            self.vert_shift_axis = np.copy(toric_nc.variables['vert_shift_axis'].data)
            self.vert_shift_wall = np.copy(toric_nc.variables['vert_shift_wall'].data)
            self.vert_triang_edge = np.copy(toric_nc.variables['vert_triang_edge'].data)
            self.vert_triang_wall = np.copy(toric_nc.variables['vert_triang_wall'].data)
            self.edge_skewdness = np.copy(toric_nc.variables['edge_skewdness'].data)
            self.dist_plafars = np.copy(toric_nc.variables['dist_plafars'].data)
            self.dist_plaant = np.copy(toric_nc.variables['dist_plaant'].data)
            self.dist_plawall = np.copy(toric_nc.variables['dist_plawall'].data)
            self.centr_elec_dens = np.copy(toric_nc.variables['centr_elec_dens'].data)
            self.centr_elec_temp = np.copy(toric_nc.variables['centr_elec_temp'].data)
            self.sep_elec_dens = np.copy(toric_nc.variables['sep_elec_dens'].data)
            self.sep_elec_temp = np.copy(toric_nc.variables['sep_elec_temp'].data)
            self.centr_ion_temp = np.copy(toric_nc.variables['centr_ion_temp'].data)
            self.sep_ion_temp = np.copy(toric_nc.variables['sep_ion_temp'].data)
            self.so_thickness = np.copy(toric_nc.variables['so_thickness'].data)
            self.so_dens_length = np.copy(toric_nc.variables['so_dens_length'].data)
            self.so_ele_temp_length = np.copy(toric_nc.variables['so_ele_temp_length'].data)
            self.so_ion_temp_len = np.copy(toric_nc.variables['so_ion_temp_len'].data)
            self.ppnei = np.copy(toric_nc.variables['ppnei'].data)
            self.ppnee = np.copy(toric_nc.variables['ppnee'].data)
            self.pptei = np.copy(toric_nc.variables['pptei'].data)
            self.pptee = np.copy(toric_nc.variables['pptee'].data)
            self.pptii = np.copy(toric_nc.variables['pptii'].data)
            self.pptie = np.copy(toric_nc.variables['pptie'].data)
            self.relef = np.copy(toric_nc.variables['relef'].data)
            self.ielef = np.copy(toric_nc.variables['ielef'].data)
            self.poynt_flux = np.copy(toric_nc.variables['poynt_flux'].data)
            self.pow_prof_elec_fw = np.copy(toric_nc.variables['pow_prof_elec_fw'].data)
            self.pow_prof_eld = np.copy(toric_nc.variables['pow_prof_eld'].data)
            self.pow_prof_ttmpe = np.copy(toric_nc.variables['pow_prof_ttmpe'].data)
            self.pow_prof_mxde = np.copy(toric_nc.variables['pow_prof_mxde'].data)
            self.pow_prof_ibwe = np.copy(toric_nc.variables['pow_prof_ibwe'].data)
            self.pow_prof_tot_elec = np.copy(toric_nc.variables['pow_prof_tot_elec'].data)
            self.pow_prof_ICfund_ions = np.copy(toric_nc.variables['pow_prof_ICfund_ions'].data)
            self.pow_prof_ICharm_ions = np.copy(toric_nc.variables['pow_prof_ICharm_ions'].data)
            self.specif_volume = np.copy(toric_nc.variables['specif_volume'].data)
            self.specif_area = np.copy(toric_nc.variables['specif_area'].data)
            self.tot_rf_current = np.copy(toric_nc.variables['tot_rf_current'].data)
            self.zeff = np.copy(toric_nc.variables['zeff'].data)
            self.prof_rf_curr = np.copy(toric_nc.variables['prof_rf_curr'].data)
            self.equil_file = np.copy(toric_nc.variables['equil_file'].data)
            self.profnt_file = np.copy(toric_nc.variables['profnt_file'].data)

        #small vars needed for stuff like power rescaling, etc.
        self.frequency = np.copy(toric_nc.variables['frequency'].data)
        self.atomic_masses = np.copy(toric_nc.variables['atomic_masses'].data)
        self.atomic_charges = np.copy(toric_nc.variables['atomic_charges'].data)
        self.concentrations = np.copy(toric_nc.variables['concentrations'].data)
        self.psi_mesh = np.copy(toric_nc.variables['psi_mesh'].data)
        self.total_power = np.copy(toric_nc.variables['total_power'].data)
        self.tot_power_eles = np.copy(toric_nc.variables['tot_power_eles'].data)
        self.fw_power_elec = np.copy(toric_nc.variables['fw_power_elec'].data)
        self.ib_power_elec = np.copy(toric_nc.variables['ib_power_elec'].data)
        self.ICfund_power_ions = np.copy(toric_nc.variables['ICfund_power_ions'].data)
        self.ICharm_power_ions = np.copy(toric_nc.variables['ICharm_pwer_ions'].data)
    
        toric_nc.close()

    #return power rescaling target values
    def pwrtarget_vals(self):

        m = self.atomic_masses
        q = self.atomic_charges
        specidx = np.arange(self.n_of_ionspec)

        #find the minority species charge mass ratio
        #first check if last species is fast alphas or minority ion
        if ((int(m[-1])==4)and(int(q[-1]==2))):
            chrgm = int(q[-2])/int(m[-2])
            minidx = specidx[-2]
        else:
            chrgm = int(q[-1])/int(m[-1])
            minidx = specidx[-1]

        #find the "pair" species with the same charge/mass ratio as the
        #minority ion.
        pairidx=0
        for i in specidx:
            spec_chrgm = int(q[i])/int(m[i])
            if ((spec_chrgm == chrgm)and(i != minidx)):
                pairindx = i

        tpw = self.total_power
        pfrac_min  = self.ICfund_power_ions[minidx]/tpw
        pfrac_pair = self.ICharm_power_ions[pairidx]/tpw
        pfrac = [pfrac_min,pfrac_pair]

        return pfrac
        

    
