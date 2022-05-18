import matplotlib.pyplot as plt
import scipy.io as sio
import imp
import h5py
import numpy as np
import os
import sys
from scipy.constants import c
from scipy.linalg import expm
from scipy.optimize import curve_fit
from scipy.misc import derivative
from itertools import product
import copy
sys.path.append("/opt/lumerical/v212/api/python") # Change this to the path on your computer
import lumapi


class mode_sol:

    def __init__(self,simulation_parameters,mode = 'TE',hide=False,close=False):

        self.simulation_parameters = simulation_parameters
        self.hide = hide
        self.close = close
        self.transverse_mode = mode
        print('Initialising MODE instance...')
        self.mode = lumapi.MODE(hide = self.hide)
        # self.mode.eval("cd('" + os.getcwd() + "');")
        self.mode.eval("cd('" + os.path.join(os.getcwd(),'LumericalScripts') + "');")
        
        print('Inputting simulation parameters...')
        self.input_parameters(self.simulation_parameters)
        print('Fitting material parameters...')
        if self.simulation_parameters['material_wg'] == "SiN (Silicon Nitride) - MuensterSpecs":
            self.mode.eval("ImportMaterials;")
        self.fit_materials()

        print('Initialisation complete')


    def close_instance(self):
        self.mode.close()

    def input_parameters(self,sim_params):
        """Input parameters into Lumerical instance"""
        for param in sim_params:

            self.mode.putv(param,sim_params[param])


    def fit_materials(self,coeffs = 3,tol = 0.01):

        """fitting effective index model for the required materials in Lumerical"""

        materials = [self.simulation_parameters['material_wg'],self.simulation_parameters['material_clad'],"Si (Silicon) - Palik","SiO2 (Glass) - Palik"]
        
        num_materials = len(materials)

        self.mode.putv('num_materials',int(num_materials))
        self.mode.putv('max_coefficients',coeffs)
        self.mode.putv('tolerance',tol)

        for i,material in enumerate(materials):
            param = f'mat{i+1}'
            self.mode.putv(param, material)
            

        self.mode.eval('fit_materials;')
        wg = self.simulation_parameters['material_wg']
        clad = self.simulation_parameters['material_clad']
        self.mode.putv('material_wg', f'{wg} Copy 1')
        self.mode.putv('material_clad', f'{clad} Copy 1')


    def find_guided(self,array,neff_array,is_TE):
        
        '''checks for guided modes either TE or TM. returns array of parameters only for guided modes
        
        Could be better to do guiding check in lumerical for frequency sweep (don't do frequency sweeps for not guided modes)

        '''

        array_checked = []

        for i,neff in enumerate(neff_array):

            if self.transverse_mode == 'TE':

                if is_TE[i] > 0.9 and  min(np.real(neff)) > 1.444:

                    array_checked.append(array[i])

            elif self.transverse_mode == 'TM':


                if is_TE[i] < 0.1 and  min(np.real(neff)) > 1.444:
                    
                    array_checked.append(array[i])

        return np.array(array_checked)

    def freq_sweep(self):
        
        """Frequency sweep performed on a single waveguide geometry"""

        print('Performing frequency sweep')        
        self.mode.eval("run_freq_sweep;")
        self.mode.eval("save;")
        print('Sweep complete')
        
        neff =np.real(self.mode.getv("neff"))
        vg = self.mode.getv("vg")
        is_TE = abs(self.mode.getv("is_TE"))
        D = self.mode.getv("D")
        f = (self.mode.getv("f")).T
        wav = c / f[0]
        neff_checked = self.find_guided(neff,neff,is_TE)
        vg_checked = self.find_guided(vg,neff,is_TE)
        D_checked = self.find_guided(D,neff,is_TE)
        
        if self.close: self.close_instance()
        
        return neff_checked,vg_checked,D_checked,wav
    def parameter_sweep(self,sim_params):
        """Sweep a design parameter eg wg width inside lumerical 
         sim_params should be dictionary with entries parameter,parameter_start,parameter_end and steps"""
        self.mode.putv("freq_sweep",0)
       
        self.input_parameters(sim_params)

        print('Performing parameter sweep')        
        self.mode.eval("parameter_sweep;")
        self.mode.eval("save;")
        print('Sweep complete')
        neff =np.real(self.mode.getv("neff"))
       
        is_TE = abs(self.mode.getv("is_TE"))
        
        if self.close: self.close_instance()
        
        neff_checked = [self.find_guided(neff[i],neff[i],is_TE[i]) for i in range(sim_params['steps'])]

        return neff_checked

        
    def parameter_freq_sweep(self,sim_params): 
        
        """ Performs frequency sweep and returns group index of only the guided modes - TE or TM specified by class attribute
        
        sim_params should be dictionary with entries parameter,parameter_start,parameter_end and steps

       """
        self.input_parameters(sim_params)

        print('Performing parameter sweep')        
        self.mode.eval("parameter_freq_sweep;")
        self.mode.eval("save;")
        print('Sweep complete')
        neff =np.real(self.mode.getv("neff"))
        vg_array = self.mode.getv("vg")
        neff_array =self.mode.getv("neff")
        D_array = self.mode.getv("D")
        is_TE = abs(self.mode.getv("is_TE"))
        f_vg = (self.mode.getv("f_vg")).T

     
        neff_checked = []
        vg_checked = []
        D_checked = []

        for i in range(sim_params['steps']):

            neff_checked.append(self.find_guided(neff_array[i],neff_array[i],is_TE[i]))
            vg_checked.append(self.find_guided(vg_array[i],neff_array[i],is_TE[i]))
            D_checked.append(self.find_guided(D_array[i],neff_array[i],is_TE[i]))
           
        
        wav = c/f_vg

        
        if self.close: self.close_instance()
            
        return neff_checked,vg_checked,D_checked,wav

    def Extract_E_field(self,freq_sweep = 0):
        """Extract Electric field profile and refractive index profile of supported modes for a single waveguide geometry"""
        print('Getting mode info...')
        self.mode.putv("freq_sweep",  freq_sweep)

        self.mode.eval("get_mode_info;")
        self.mode.eval("save;")
        vg_checked = []

        neff =np.real(self.mode.getv("neff"))
        # print(neff)
        is_TE = self.mode.getv("is_TE")
        if freq_sweep:
            vg_array = self.mode.getv("vg")
            vg_checked = self.find_guided(vg_array,neff,is_TE)
        x_array = self.mode.getv("x_array")[:,0]
        y_array = self.mode.getv("y_array")[:,0]

        # print(x_array.shape)
        if self.transverse_mode == 'TE':

            self.E_field =self.mode.getv("Ex")

        elif self.transverse_mode == 'TM':

            self.E_field = self.mode.getv("Ey")

        index_profile_temp =self.mode.getv("index_profile")
        index_profile = np.array(index_profile_temp[:,:,0,0]).T
        # print(E_field.shape)
        # print(neff)
        E_filtered = self.find_guided(self.E_field[:,:,:,0,0],neff,is_TE)
        # print(E_filtered.shape)
        array_dimension = []
        array_dimension.append(index_profile.shape[0])
        array_dimension.append(index_profile.shape[1])

        if self.close: self.close_instance()
            
        return E_filtered ,index_profile,x_array,y_array,vg_checked

