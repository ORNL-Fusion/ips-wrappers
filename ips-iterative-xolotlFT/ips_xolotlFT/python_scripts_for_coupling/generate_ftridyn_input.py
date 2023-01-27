"""
Provides classes and functions to generate F-TRIDYN input files.
"""

from __future__ import print_function
import numpy as np
#import matplotlib.pyplot as plt
import os

class species:
	"""
	Stores species-relevant parameters for F-TRIDYN input.

	Args:
		ZZ (int): Atomic number.
		M (float): Mass in amu.
		BE (float): Bulk binding energy in eV.
		ED (float): Displacement energy in eV.
		EF (float): Threshold (cutoff) energy in eV.
		QU (float): Proportion of species in target.
		DNS0 (float): Number density of species in inverse Angstroms cubed.
		CK (float): Electronic stopping correction factor.
		E0 (float): Incident energy in eV.
		ALPHA0 (float): Incident energy in degrees.
		QUBEAM (float): Proportion of species in beam.
		QUMAX (float): Max proportion of species allowed in target.
		SBV (float): Surface binding energy in eV.
	"""
	def __init__(self, ZZ, M, BE, ED, EF, QU, DNS0, CK, E0, ALPHA0, QUBEAM, QUMAX, SBV):
		self.ZZ = ZZ					#Atomic number/charge
		self.M = M						#Mass in [u]
		self.BE = BE					#Binding energy in [eV]
		self.ED = ED					#Displacement energy in [eV]
		self.EF = EF					#Cutoff energy in [eV]
		self.QU = QU					#Concentration of species in homogeneous target
		self.DNS0 = DNS0				#Number density in [at/A^3]
		self.CK = CK					#Electronic stopping correction factor (0.0<CK<1.0)
		self.E0 = E0					#Incident energy in [eV]
		self.ALPHA0 = ALPHA0			#Incident angle in [deg] (0.0<ALPHA0<90.0)
		self.QUBEAM = QUBEAM			#Concentration in beam (0.0<QUBEAM<1.0)
		self.QUMAX = QUMAX				#Maximum allowed concentration in target (0.0<QUMAX<1.0)
		self.SBV = SBV					#Surface binding energy in [eV]
	#end def __init__
#end class species

class species_lookup_table:
	"""
	Provides functions to lookup and generate species objects for F-TRIDYN
	input. Searches for file 'table1.txt' and stores the first section.
	"""
	def __init__(self):
		cwd = os.path.dirname(os.path.realpath(__file__))
		self.table1 = open(cwd+'/'+'table1.txt','r').readlines()[11:130]
		self.species_name = []
		self.ZZ = []
		self.M = []
		self.BE = []
		self.ED = []
		self.EF = []
		self.DNS0 = []
		self.CK = []
		self.SBV = []
		for i in range(len(self.table1)):
			line = self.table1[i]
			self.species_name.append(line[0:5].strip())
			self.ZZ.append(float(line[6:10].strip()))
			self.M.append(float(line[11:22].strip()))
			self.BE.append(float(line[114:120].strip()))
			self.ED.append(float(line[53:59].strip()))
			self.EF.append(float(line[62:68].strip()))
			self.DNS0.append(float(line[32:41].strip()))
			self.SBV.append(float(line[43:49].strip()))
		#end for
	#end def __init__

	def find_species(self,name,QU=0.0,CK=1.0,E0=0.0,ALPHA0=0.0,QUBEAM=0.0,QUMAX=0.0):
		"""
		Lookup species parameters from included data file and return an appropriate
		species object for F-TRIDYN input.

		Args:
			name (:obj:'str'): species name (i.e., chemical symbol).
			QU (float, optional): portion of species in target.
			CK (float, optional): Electronic stopping correction factor.
			E0 (float, optional): Incident energy of species in eV.
			ALPHA0 (float, optional): Incident angle of species in degrees.
			QUBEAM (float, optional): Portion of species in beam.
			QUMAX (float, optional): Maximum allowed proportion in target.
		"""
		species_index = self.species_name.index(name)
		return species(ZZ=self.ZZ[species_index],
			M=self.M[species_index], BE=self.BE[species_index],
			ED=self.ED[species_index], EF=self.EF[species_index],
			QU=QU, DNS0=self.DNS0[species_index], CK=CK, E0=E0, ALPHA0=ALPHA0,
			QUBEAM=QUBEAM, QUMAX=QUMAX, SBV=self.SBV[species_index])
	#end def find_species
#end def species_lookup

class sim_params:
	"""
	Holds simulation parameters for F-TRIDYN simulations and provides methods for
	generating and producing input files.

	Args:
		name (:obj:'str'): a 4-character name that represents the simulation (e.g., 'He_W')
		sim_num (int): Number of simulation. Useful for ordering multiple simulations with same name.
		IFOUT (int): Number of projectiles after which to write info to terminal.
		NH (int): Number of projectile histories (trajectories) to perform.
		IDOUT (int): Number of projectiles after which to produce integral output.
		IQOUT (int): Number of projectiles after which to produce profile files (.PR##).
		NCP (int): Number of speicies in simulation.
		IDREL (int): When > 0, relaxation is suppressed. When < 0, relaxation and cascades are suppressed.
		IQ0 (int): When == 0, homogeneous initial composition. When != 0, inhomogeneous initial composition from .LAY file.
		IRC0 (int): When < 0, subthreshold atoms are free. When >= 0, subthreshold atoms are bound.
		IRAND (int): Seed for PRNG.
		JSP1 (int): Suppression of recoils of species JSP1 to JSP2.
		JSP2 (int): Suppression of recoils of species JSP1 to JSP2.
		JFRP (int): Frenkel Pair generation of species 1 to JFRP.
		JNRM (int): Normalization of print output.
		FLC (float): Fluence for dynamic composition simulation.
		INEL (int): Electronic stopping model. 1: inelastic nonlocal; 2: local; 3: equipartition
		IWC (int): Order of weak collisions (i.e., how many possible weak collisions per collision event)
		IDIFF (int): Extra particle handling. 1: Reemission; 2: "Diffusion"
		CSBE (int): Surface binding energy model.
		ANGMN (float): Minimum angle to include in particle list output.
		ANGMX (float): Maximum angle to include in particle list output.
		TT (float): Depth of surface.
		TTDYN (float): Depth of dynamic surface. Defaults to 20% more than TT.
		NQX (int): Number of composition layers.
		DSF (float): Averaging depth for surface composition output.
		IQXN (int): Limiting depth intervals for profile output.
		IQXX (int): Unused.
		IMCP (int): Species for which moments are calculated.
		surf_name (:obj:'str'): Name of surface file.
		species_list (:obj:'list'): List of species in simulation.
		output_pka (int, optional): Whether to output PKA list. Defaults to 0. [0,1]
		output_ska (int, optional): Whether to output SKA list. Defaults to 0. [0,1]
		output_prj (int, optional): Whether to output implanted projectile list. Defaults to 1. [0,1]
		output_out (int, optional): Whether to output human-readable output. Defaults to 1. [0,1]
		output_srfc (int, optional): Whether to output surface composition. Defaults to 0. [0,1]
		output_spyl (int, optional): Whether to output sputtering yield. Defaults to 0. [0,1]
		output_ardn (int, optional): Unused output. Defaults to 0. [0,1]
		output_srrs (int, optional): Unused output. Defaults to 0. [0,1]
		output_reem (int, optional): Whether to output reemitted particles. Defaults to 0. [0,1]
		output_splst (int, optional): Whether to output sputtered list. Defaults to 1. [0,1]
		output_rflst (int, optional): Whether to output reflected list. Defaults to 0. [0,1]
		output_trlst (int, optional): Whether to output transmitted list. Defaults to 0. [0,1]
		output_simple (int, optional): Whether to output simple output. Defaults to 0. [0,1]
		output_fp (int, optional): Whether to output Frenkel Pair list. Defaults to 0. [0,1]
		output_prof (int, optional): Whether to output composition profiles. Defaults to 0. [0,1]

	"""
	def __init__(self,name,sim_num,IFOUT,NH,IDOUT,IQOUT,NCP,IDREL,IQ0,IRC0,IRAND,
		JSP1,JSP2,JFRP,JNRM,FLC,INEL,IWC,IDIFF,CSBE,ANGMN,ANGMX,TT,TTDYN,NQX,DSF,
		IQXN,IQXX,IMCP,surf_name,species_list, output_pka=0, output_ska=0, output_prj=1,
		output_out=1,output_srfc=0,output_spyl=0,output_ardn=0,output_srrs=0,output_reem=0,
		output_splst=1,output_rflst=0,output_trlst=0,output_simple=0,output_fp=0,output_prof=0):
		self.name = name 						#4-character file name, e.g., He_W, C_C_, NaNa, etc.
		self.sim_num = str(sim_num).zfill(4) 	#4-number simulation number, e.g., 0001, 0002
		self.IFOUT = int(IFOUT) 				#Number of projectiles after which info is written to terminal
		self.NH = int(NH)						#Number of incident histories (i.e. projectiles)
		self.IDOUT = int(IDOUT) 				#Integral data ouput after every IDOUT'th history
		self.IQOUT = int(IQOUT) 				#Additional .PR## profile file output after every IQOUT'th history
		self.NCP = int(NCP) 					#Number of components (i.e. species)
		self.IDREL = int(IDREL) 				#IDREL>0:Suppression of relaxation, IDREL<0:Suppression of relaxation and cascades
		self.IQ0 = int(IQ0) 					#IQ0=0:homogeneous initial composition, IQ0!=0:inhomogeneous initial composition from .LAY
		self.IRC0 = int(IRC0) 					#IRC0<0: Subthreshold atoms free, IRC0>=0: Subthreshold atoms bound
		self.IRAND = int(IRAND) 				#Random seed
		self.JSP1 = int(JSP1) 					#Suppression of recoils of type JSP1,...,JSP2
		self.JSP2 = int(JSP2) 					#JSP1=0: JSP2 ignored, all recoils traced
		self.JFRP = int(JFRP) 					#Frenkel Pair generation for JFRP,...,NCP
		self.JNRM = int(JNRM) 					#Normalzation of print output
		self.FLC = FLC 							#Fluence [1e16/cm2]
		self.INEL = int(INEL) 					#INEL=[1,2,3]:[inelastic nonlocal, local, equipartition]
		self.IWC = int(IWC) 					#Max order of weak projectile-target collisions {1,2,3}
		self.IDIFF = int(IDIFF) 				#IDIFF=[0,1]:[reemission,"diffusion"]
		self.CSBE = CSBE 						#Surface binding energy model (see documentation)
		self.ANGMN = ANGMN 						#Minimum angle for sputtered energy/angle distributions
		self.ANGMX = ANGMX 						#Maximum angle for sputtered energy/angle distributions
		self.TT = TT 							#Depth of surface
		self.TTDYN = TTDYN 						#Depth of dynamic surface (>=TT)
		self.NQX = NQX 							#Number of composition layers within F-TRIDYN
		self.DSF = DSF 							#Averaging depth for surface composition
		self.IQXN = IQXN 						#Limiting depth intervals for .PR## profile output
		self.IQXX = IQXX						#
		self.IMCP = int(IMCP) 					#Component for which moments should be calcluated

		self.surf_name = surf_name 				#Name of surface file
		self.species_list = species_list 		#List of species objects

		self.output_pka = int(output_pka) 		#PKA output from F-TRIDYN
		self.output_ska = int(output_ska) 		#SKA output from F-TRIDYN
		self.output_prj = int(output_prj) 		#Projectile output from F-TRIDYN

		self.output_out = int(output_out) 		#Human-readable output file
		self.output_srfc = int(output_srfc) 	#
		self.output_spyl = int(output_spyl) 	#Sputtering yield calculation
		self.output_ardn = int(output_ardn) 	#
		self.output_srrs = int(output_srrs) 	#
		self.output_reem = int(output_reem) 	# Re-emission
		self.output_splst = int(output_splst) 	#Sputtered particle list
		self.output_rflst = int(output_rflst) 	#Reflected particle list
		self.output_trlst = int(output_trlst) 	#Transmitted particle list
		self.output_simple = int(output_simple) #Simple output file
		self.output_fp = int(output_fp) 		#Frenkel Pair creation site list
		self.output_prof = int(output_prof) 	#Dynamic composition profiles
		#todo: error checking, NCP=length(species_list), all parameters within bounds, etc.
	#end def __init__

	def print_input_file(self):
		"""
		Generates an F-TRIDYN input file with a .IN extension.

		Input file format is as follows:
		NAMExxxx.IN
		IFOUT
		NH  IDOUT  IQOUT  NCP  IDREL  IQ0  IRC0  IRAND  JCP1  JCP2  JFRP  JNRM
		FLC  INEL  IWC  IDIFF  CSBE  ANGMN  ANGMX
		TT  TTDYN  NQX  DSF  IQXN  IQXX  IMCP
		surfname.txt output_flags ...
		ZZ(1)  M(1)  BE(1)  ED(1)  EF(1)  QU(1)  DNS0(1)  CK(1)
		E0(1)  ALPHA0(1)  QUBEAM(1)  QUMAX(1)  SBV(1)  SBV(1)
		ZZ(2)  M(2)  BE(2)  ED(2)  EF(2)  QU(2)  DNS0(2)  CK(2)
		E0(2)  ALPHA0(2)  QUBEAM(2)  QUMAX(2)  SBV(2)  SBV(2)
		ZZ(3)...
		ZZ(NCP)...
		"""
		input_file = open(self.name+self.sim_num+'.IN','w+')
		print(self.name+self.sim_num+'.IN', file=input_file)
		print(self.IFOUT,file=input_file)
		print(self.NH, self.IDOUT, self.IQOUT, self.NCP, self.IDREL, self.IQ0,
			self.IRC0, self.IRAND, self.JSP1, self.JSP2, self.JFRP, self.JNRM,
			file=input_file)
		print(self.FLC, self.INEL, self.IWC, self.IDIFF, self.CSBE, self.ANGMN,
			self.ANGMX, file=input_file)
		print(self.TT, self.TTDYN, self.NQX, self.DSF, self.IQXN, self.IQXX,
			self.IMCP, file=input_file)
		print(self.surf_name,self.output_pka, self.output_ska, self.output_prj,
			self.output_out,self.output_srfc,self.output_spyl,self.output_ardn,
			self.output_srrs,self.output_reem,self.output_splst,self.output_rflst,
			self.output_trlst,self.output_simple,self.output_fp,self.output_prof,
			file=input_file)
		for i in range(self.NCP):
				current_species = self.species_list[i]

				print(current_species.ZZ, current_species.M, current_species.BE,
					current_species.ED, current_species.EF, current_species.QU,
					current_species.DNS0, current_species.CK, file=input_file)

				print(current_species.E0, current_species.ALPHA0,
					current_species.QUBEAM, current_species.QUMAX, end='',
					file=input_file)

				for j in range(self.NCP):
					print(' '+str(current_species.SBV), end='', file=input_file)
				#end for

				print(file=input_file)
		#end for
	#end def print_input_file
#end class simulation

class fractal_surface:
	"""
	Stores and provides methods for generating Mandelbrot-style fractal surfaces.

	Args:
		FD (float, optional): fractal dimension of generated fractal surface.
			Defaults to 1.0, a smooth surface.
		width (float, optional): Total width of fractal surface in Angstroms.
			Defaults to 400 Angstroms.
		iterations (int, optional): Number of generations to iterate. Defaults to 1.

	"""
	def __init__(self, FD=1.0, width=400.0, iterations=1):
		self.FD = FD
		self.width = width
		if FD==1.0:
			self.iterations = 1
		else:
			self.iterations = iterations
		self.fgen()
		self.name = str(1)+'p'+str(int((FD-0.99999)*1000.0))
	#end def __init__

	def fgen(self,shape=[0,1,0,-1,0,-1,0,1,0]):
		"""
		Given a shape, produce a generated fractal using Mandelbrot's fractal
		generator method, whereby an initial shape is transformed onto each
		segment of itself over a specified number of generations to produce a
		curve with an arbitrary fractal dimension given an initial shape type.

		Args:
			shape (:obj:'list', optional) This function defines the initial
			fractal generator. Values in the list are 0 (horizontal line
			segment), 1 (upward angled line segment) and -1 (downward angled
			line segment). Sum of shape should be zero for consistency. Angles
			are calculated from FD.

			Defaults to [0,1,0,-1,0,-1,0,1,0]. This produces surfaces that look
			like the following illustration:

			y=0       y=width
			|  _      |
			|_/ \_   _|_ x=0
			|     \_/ |

			As another example of a different shape, the shape [0,1,0,-1,0]
			looks like the following illustration:

			y=0     y=width
			|       |
			|_/\_  _|_ x=0
			|    \/ |
		"""
		b = np.sum(np.abs(shape))
		a = np.size(shape)-b
		beta = np.arccos((-a+(a+b)**(1.0/self.FD))/b)
		L = 1.0/(a+b*np.cos(beta))
		num_gen_seg = np.size(shape)
		self.num_gen_points = num_gen_seg+1

		self.x = np.zeros(self.num_gen_points)
		self.y = np.zeros(self.num_gen_points)
		segX = np.zeros(self.num_gen_points)
		segY = np.zeros(self.num_gen_points)
		self.gx = np.zeros(self.num_gen_points)
		self.gy = np.zeros(self.num_gen_points)

		self.x[0] = 0.0
		self.y[0] = 0.0

		for i in range(1,self.num_gen_points):
			self.x[i] = self.x[i-1] + L * np.cos(beta*shape[i-1])
			self.y[i] = self.y[i-1] + L * np.sin(beta*shape[i-1])
		#end for

		self.gx = self.x
		self.gy = self.y

		for i in range(2,self.iterations+1):
			counter = -1
			storeX = self.x
			storeY = self.y
			self.x = np.zeros(num_gen_seg**(i)+1)
			self.y = np.zeros(num_gen_seg**(i)+1)

			for j in range(np.size(storeX)-1):
				segX[0] = storeX[j]
				segY[0] = storeY[j]
				gamma = np.arctan2((storeY[j+1]-segY[0]),(storeX[j+1]-segX[0]))

				for k in range(1,num_gen_seg+1):
					coordX = L**i*np.cos(beta*shape[k-1])
					coordY = L**i*np.sin(beta*shape[k-1])
					segX[k] = segX[k-1] + coordX*np.cos(gamma) - coordY*np.sin(gamma)
					segY[k] = segY[k-1] + coordX*np.sin(gamma) + coordY*np.cos(gamma)
				#end for k

				if(j<np.size(storeX)-2):
					segTrim = -1
				else:
					segTrim = 0
				#end if

				for n in range(0,np.size(segX)+segTrim):
					counter = counter + 1
					self.x[counter] = segX[n]
					self.y[counter] = segY[n]
				#end for n
			#end for j
		#end for i
		self.num_surf_points = np.size(self.x)
		self.x,self.y = self.width*self.x,self.width*self.y #scale the final surface
	#end def fgen

	def print_fractal_surface_file(self):
		"""
		Print fractal surface file for F-TRIDYN input.

		File format is as follows:
		<number of generator points> <number of iterations> <number of points>
		<fractal dimension>
		<generator x points> <generator y points>
		...
		<fractal x points> <fractal y points>
		...

		"""
		surface_file = open(self.name+'.surf','w+')

		print(self.num_gen_points, self.iterations, self.num_surf_points, file=surface_file)
		print(self.FD,file=surface_file)

		for i in range(self.num_gen_points):
			print(self.gx[i], self.gy[i], file=surface_file)
		#end for

		for i in range(self.num_surf_points):
			print(self.x[i], self.y[i], file=surface_file)
		#end for

	#end def print_fractal_surface_file
#end class fractal_surface

def He_W_xolotl(sim_number=1, number_histories=1000000,
	incident_energy=200.0, incident_angle=0.0, fractal_dimension=1.0,
	IQ0=0, number_layers=100, depth=200.0):

	#Lookup species physical parameters from table
	lookup_table = species_lookup_table()
	Helium = lookup_table.find_species('He', E0=incident_energy, ALPHA0=incident_angle, QUBEAM=1.0, QUMAX=0.0)
	Tungsten = lookup_table.find_species('W', QU=1.0, QUMAX=1.0)

	#Collect species into species_list
	species_list = [Helium, Tungsten]
	number_species = np.size(species_list)

	#Create fractal and print surface file to 1p###.surf
	surface = fractal_surface(FD=fractal_dimension, width=400.0, iterations=3)
	surface.print_fractal_surface_file()

	#Define simulation parameters for F-TRIDYN
	simulation_parameters = sim_params(name='He_W', sim_num=sim_number, IFOUT=int(number_histories/20),
		NH=int(number_histories), IDOUT=int(number_histories/5), IQOUT=int(number_histories/5), NCP=number_species, IDREL=1,
		IQ0=IQ0, IRC0=-1, IRAND=12855897, JSP1=0, JSP2=1, JFRP=1, JNRM=1, FLC=1E-16, INEL=1,
		IWC=3, IDIFF=1, CSBE=0, ANGMN=0, ANGMX=90, TT=depth, TTDYN=1.2*depth, NQX=number_layers,
		DSF=100.0, IQXN=0, IQXX=250, IMCP=0, surf_name=surface.name+'.surf', species_list=species_list,
		output_pka=0, output_ska=0, output_prj=1)

	#Print simulation parameters to <name><sim_number>.IN
	simulation_parameters.print_input_file()
#end def He_W

def Prj_Tg_xolotl(sim_number=1, number_histories=1000000,
                  incident_energy=200.0, incident_angle=0.0, fractal_dimension=1.0,
                  IQ0=0, number_layers=100, depth=200.0, projectile_name='He',target1_name='W', target2_name='',target3_name='',target4_name='',**kwargs):

        #Lookup species physical parameters from table
        lookup_table = species_lookup_table()
        projectile=lookup_table.find_species(projectile_name, E0=incident_energy, ALPHA0=incident_angle, QUBEAM=1.0, QUMAX=0.0)
        target1=lookup_table.find_species(target1_name, QU=1.0, QUMAX=1.0)

        #Collect species into species_list
        species_list = [projectile,target1]

		#Update species list with additional parameters
        for prop in ["ED_He", "EF_He", "E0_He", "ALPHA0_He", "ED_W", "EF_W", "SBV_W"]:
            if prop in kwargs.keys():
                n = prop.endswith("_W")
                setattr(species_list[n], prop, kwargs[prop])

        #add target species if other than empty:
        #currently max of 5 species handled by F-Tridyn
        additional_targets=[target2_name, target3_name, target4_name]

        for  sp in additional_targets:
                if sp:
                        add_target=lookup_table.find_species(sp, QU=1.0, QUMAX=1.0)
                        species_list.append(add_target)


        number_species = np.size(species_list)

        #Create fractal and print surface file to 1p###.surf
        surface = fractal_surface(FD=fractal_dimension, width=400.0, iterations=3)
        surface.print_fractal_surface_file()

        #string of 4-characters required for name
        if (len(projectile_name)+len(target1_name))==4:
                parameterNameString=projectile_name+target1_name #already 4 char, e.g. HeTa
        elif (len(projectile_name)+len(target1_name))==3:
                parameterNameString=projectile_name+'_'+target1_name #e.g.: He_W
        else: #len(parameterNameString) ==2
                parameterNameString=projectile_name+'_'+target1_name+'_' #e.g., W_W -> W_W_

        #Define simulation parameters for F-TRIDYN
        simulation_parameters = sim_params(name=parameterNameString, sim_num=sim_number, IFOUT=int(number_histories/20),
                NH=int(number_histories), IDOUT=int(number_histories/5), IQOUT=int(number_histories/5), NCP=number_species, IDREL=1,
                IQ0=IQ0, IRC0=-1, IRAND=12855897, JSP1=0, JSP2=1, JFRP=1, JNRM=1, FLC=1E-16, INEL=1,
                IWC=3, IDIFF=1, CSBE=0, ANGMN=0, ANGMX=90, TT=depth, TTDYN=1.2*depth, NQX=number_layers,
                DSF=100.0, IQXN=0, IQXX=250, IMCP=0, surf_name=surface.name+'.surf', species_list=species_list,
                output_pka=0, output_ska=0, output_prj=1)

        #Print simulation parameters to <name><sim_number>.IN
        simulation_parameters.print_input_file()
#end def Prj_Tg

def beam_and_target(name,beam_species,target_species,sim_number=1,
	number_histories=1E5,incident_energy=100.0,incident_angle=0.0,
	fractal_dimension=1.0,width=200.0,depth=200.0):
	"""
		A function that generates an input file given a monoenergetic beam on a
		homogeneous target.

		Args:
			name (:obj:'str'): a 4-character name that represents the simulation (e.g., 'He_W')
			beam_species (:obj:'str'): Chemical symbol of beam species.
			target_species (:obj:'str'): Chemical symbol of target species.
			sim_number (int, optional): Number of simulation. Useful for ordering multiple simulations with same name.
			number_histories (int, optional): Number of histories (trajectories). Defaults to 1E5.
			incident_energy (float, optional): Incident beam energy. Defaults to 100.0. [eV]
			incident_angle (float, optional): Incident angle of beam. Defaults to 0.0. [degrees]
			fractal_dimension (float, optional): FD of surface. Defaults to 1.0 (smooth).
			width (float, optional): Total width of surface. [A]
			depth (float, optional): Total depth of surface. [A]

		Returns:
			simulation_parameters (:obj:'sim_params'): sim_params object.

		Prints:
			A single F-TRIDYN input file (extension: .IN) for use in F-TRIDYN simulations.

	"""

	lookup_table = species_lookup_table()
	beam = lookup_table.find_species(beam_species, E0=incident_energy, ALPHA0=incident_angle, QUBEAM=1.0, QUMAX=0.0)
	target = lookup_table.find_species(target_species, QU=1.0, QUMAX=1.0)
	species_list = [beam,target]
	number_species = np.size(species_list)

	surface = fractal_surface(FD=fractal_dimension, width=width, iterations=3)
	surface.print_fractal_surface_file()

	simulation_parameters = sim_params(name=name, sim_num=sim_number, IFOUT=number_histories/20,
		NH=number_histories, IDOUT=number_histories/5, IQOUT=number_histories/5, NCP=number_species, IDREL=0,
		IQ0=0, IRC0=-1, IRAND=12855897, JSP1=0, JSP2=1, JFRP=1, JNRM=1, FLC=1e-16, INEL=1,
		IWC=3, IDIFF=0, CSBE=0, ANGMN=0, ANGMX=90, TT=depth, TTDYN=1.2*depth, NQX=500,
		DSF=100.0, IQXN=0, IQXX=250, IMCP=0, surf_name=surface.name+'.surf', species_list=species_list,
		output_pka=0, output_ska=0, output_prj=1)
	simulation_parameters.print_input_file()

	return simulation_parameters
#end def beam_and_target

def iterate_SBV(name, beam_species, target_species, incident_energy, incident_angle, depth, SBV_number, SBV_change):
	number_histories = 1E5
	lookup_table = species_lookup_table()
	beam = lookup_table.find_species(beam_species, E0=incident_energy, ALPHA0=incident_angle, QUBEAM=1.0, QUMAX=0.0)
	target = lookup_table.find_species(target_species, QU=1.0, QUMAX=1.0)
	SBV_values = target.SBV * np.arange(1.0-SBV_change,1.0+SBV_change,(2.0*SBV_change / SBV_number))
	print(SBV_values)
	surface = fractal_surface(FD=1.0, width=200, iterations=1)
	surface.print_fractal_surface_file()

	for i in range(SBV_number):
		target.SBV = SBV_values[i]
		species_list = [beam,target]
		number_species = np.size(species_list)
		simulation_parameters = sim_params(name=name, sim_num=i+1, IFOUT=number_histories/20,
			NH=number_histories, IDOUT=number_histories/5, IQOUT=number_histories/5, NCP=number_species, IDREL=0,
			IQ0=0, IRC0=-1, IRAND=12855897, JSP1=0, JSP2=1, JFRP=1, JNRM=1, FLC=1.0E2, INEL=1,
			IWC=3, IDIFF=0, CSBE=0, ANGMN=0, ANGMX=90, TT=depth, TTDYN=1.2*depth, NQX=500,
			DSF=100.0, IQXN=0, IQXX=250, IMCP=0, surf_name=surface.name+'.surf', species_list=species_list,
			output_pka=0, output_ska=0, output_prj=1)
		simulation_parameters.print_input_file()

if __name__ == '__main__':
	#How to iterate over parameters
	#for i in range(3):
	#	He_W_xolotl(sim_number='000'+str(i+1), fractal_dimension=1.0+0.1*i,number_histories=10000)
	#end for
	#beam_species = 'He'
	#target_species = 'W'
	#energy = -1
	#beam_and_target('He_W',beam_species, target_species, sim_number=1,
	#	number_histories=1E6, incident_energy=energy,depth=0.5E5,incident_angle=76.0)
	#end for
	beam_and_target('He_W','He','W')
