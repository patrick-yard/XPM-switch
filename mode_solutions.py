
import numpy as np
import os
import sys
from scipy.constants import c
from itertools import product
version = 'v212'
sys.path.append(f"C:\\Program Files\\Lumerical\\{version}\\api\\python")
import lumapi

class Solver:

    def __init__(self, simulation_parameters: dict, mode='TE', hide=False, close=False):
        '''Class to interface python and Lumerical MODE (EME) functionality. 

        Allows Lumerical simulations to be performed with python specified parameters.

        Allows eigenmode of straight waveguides to be determined and to perform frequency and parameter sweeps 

        Also can simulate eigenmodes of two coupled waveguides - should expand this to N waveguides'''

        self.simulation_parameters = simulation_parameters
        self.hide = hide
        self.close = close
        self.transverse_mode = mode
        self.initialise_mode()

    def initialise_mode(self):
        '''Initialise an instance of Lumerical mode solutions, 

            input all current simulations paramters, 

            import material if necessary - i.e. if SiN is to be used

            re-fit dispersion curves of core and cladding materials'''

        print('Initialising MODE instance...')
        self.mode = lumapi.MODE(hide=self.hide)
        self.mode.eval(
            "cd('" + os.path.join(os.getcwd(), "LumericalFiles") + "');")
        print('Inputting simulation parameters...')
        self.input_parameters(self.simulation_parameters)
        print('Fitting material parameters...')
        if self.simulation_parameters['material_wg'] == "SiN (Silicon Nitride) - MuensterSpecs":
            self.mode.eval("ImportMaterials;")
        self.fit_materials()

        print('Initialisation complete')

    def close_instance(self):
        '''Close current mode instance'''

        self.mode.close()

    def input_parameters(self, sim_params: dict):
        '''Input dictionary of simulation parameters into current mode instance'''

        for param in sim_params:

            self.mode.putv(param, sim_params[param])

    def fit_materials(self, coeffs=3, tol=0.01):
        '''Fit dispersion curves of waveguide, cladding and substrate materials'''

        materials = [self.simulation_parameters['material_wg'],
                     self.simulation_parameters['material_clad'], "Si (Silicon) - Palik", "SiO2 (Glass) - Palik"]

        num_materials = len(materials)

        self.mode.putv('num_materials', int(num_materials))
        self.mode.putv('max_coefficients', coeffs)
        self.mode.putv('tolerance', tol)

        for i, material in enumerate(materials):
            param = f'mat{i+1}'
            self.mode.putv(param, material)

        self.mode.eval('fit_materials;')
        wg = self.simulation_parameters['material_wg']
        clad = self.simulation_parameters['material_clad']
        self.mode.putv('material_wg', f'{wg} Copy 1')
        self.mode.putv('material_clad', f'{clad} Copy 1')

    def read_params(self, params: dict) -> dict:
        '''Read desired parameters from current MODE instance

        input blank dictionary with keys corresponding to the desired parameters

        returns dictionary filled with simulated results '''

        for param in params.keys():

            params[param] = self.mode.getv(param)

        return params

    def find_guided(self, array: list, neff_array: list, is_TE: list) -> np.array:
        '''checks for guided modes either TE or TM.

        Inputs desired array along with matching arrray of effective index and TE fraction values.
        Returns array of parameters only for guided modes

        Could be better to do guiding check in lumerical for frequency sweep (don't do frequency sweeps for not guided modes)

        '''

        array_checked = []

        for i, neff in enumerate(neff_array):

            if self.transverse_mode == 'TE':

                if is_TE[i] > 0.9 and min(np.real(neff)) > 1.444:

                    array_checked.append(array[i])

            elif self.transverse_mode == 'TM':

                if is_TE[i] < 0.1 and min(np.real(neff)) > 1.444:

                    array_checked.append(array[i])

        return np.array(array_checked)

    def check_many(self, params: dict, freq_sweep=False) -> dict:
        '''Find_guided for many parameters


        if a parameter sweep then neff_array and is_TE should be multidimensional containing the values for each sweep point

        if no parameter sweep they need only be a single list'''

        neff_array = params['neff']
        is_TE = params['is_TE']

        checked_params = dict()

        for key, array in params.items():
    
            if len(array.shape) == 2 or freq_sweep:

                checked_arrays_temp = self.find_guided(array, neff_array, is_TE)

            else:
                checked_arrays_temp = []

                for i, arr in enumerate(array):

                    checked_arrays_temp.append(
                        self.find_guided(arr, neff_array[i], is_TE[i]))

            checked_params[key] = np.array(checked_arrays_temp)

        return checked_params

    def freq_sweep(self) -> dict:
        '''Run a frequency sweep for a single waveguide

        Returns simulated parameters for only the guided modes 

       TODO: generalise for N-mode random walks '''

        print('Performing frequency sweep')
        self.mode.eval("run_freq_sweep;")
        self.mode.eval("save;")
        print('Sweep complete')

        res_dict = {'neff': None,
                    'vg': None,
                    'is_TE': None,
                    'D': None,
                    'f': None}

        res_dict = self.read_params(res_dict)

        f = res_dict['f'].T[0]
        res_dict.pop('f')

        checked_params = self.check_many(res_dict, freq_sweep=True)
        checked_params['wav'] = c / f
        checked_params['f'] = f

        if self.close:
            self.close_instance()

        return checked_params

    def parameter_sweep(self, sim_params: dict) -> dict:
        '''Sweep a chosen parameter for a single waveguide in Lumerical. 

        input sim_params should be dictionary with entries parameter,parameter_start,parameter_end and steps.'''

        self.mode.putv("freq_sweep", 0)

        self.input_parameters(sim_params)

        print('Performing parameter sweep')
        self.mode.eval("LumericalFiles\\parameter_sweep;")
        self.mode.eval("save;")
        print('Sweep complete')

        res_dict = {'neff': None,
                    'is_TE': None}

        res_dict = self.read_params(res_dict)

        checked_params = self.check_many(res_dict)

        return checked_params

    def parameter_freq_sweep(self, sim_params: dict) -> dict:
        """ Performs parameter sweep in lumercial with a frequency sweep at each point.

        sim_params should be dictionary with entries parameter,parameter_start,parameter_end and steps

       """
        self.input_parameters(sim_params)

        print('Performing parameter sweep')
        self.mode.eval("parameter_freq_sweep;")
        self.mode.eval("save;")
        print('Sweep complete')

        res_dict = {'neff': None,
                    'vg': None,
                    'is_TE': None,
                    'D': None,
                    'f_vg': None}

        res_dict = self.read_params(res_dict)

        res_dict['f_vg'] = res_dict['f_vg'].T

        checked_params = self.check_many(res_dict)

        checked_params['wav'] = c / res_dict['f'][0]
        checked_params['f'] = res_dict['f'][0]
        if self.close:
            self.close_instance()

        return checked_params

    def Extract_E_field(self, freq_sweep=0):
        '''Simulate E field distribution for  a single waveguide

        Input: freq_sweep: binary value which decides whether a frequency sweep is performed

        Returns: dictionary containing E field of all guided modes, 

        index profile of the waveguide structure along with x and y dimensions'''

        print('Getting mode info...')
        self.mode.putv("freq_sweep",  freq_sweep)
        self.mode.eval("get_mode_info;")
        self.mode.eval("save;")
      
        res_dict = {'neff': None,
                    'is_TE': None}
        if freq_sweep:
            res_dict.update({'vg': None})

        res_dict = self.read_params(res_dict)

        field_dict = {'x_array': None,
                      'y_array': None,
                      'index_profile': None}

        if self.transverse_mode == 'TE':
            field_str = 'Ex'
            field_dict.update({field_str: None})

        elif self.transverse_mode == 'TM':
            field_str = 'Ey'
            field_dict.update({field_str: None})

        field_dict = self.read_params(field_dict)

        checked_params = self.check_many(res_dict)

        checked_params['x_array'] = field_dict['x_array'][:, 0]
        checked_params['y_array'] = field_dict['y_array'][:, 0]        
        checked_params['index_profile'] = np.array( field_dict['index_profile'][:, :, 0, 0]).T

        checked_params['E_field'] = self.find_guided(field_dict[field_str][:, :, :, 0, 0], res_dict['neff'], res_dict['is_TE'])
        
        if self.close:
            self.close_instance()

        return checked_params

    # def parameter_sweep_DC(self, sim_params):

    #     self.input_parameters(sim_params)

    #     print('Performing parameter sweep')
    #     self.mode.eval("DC_modesolver;")
    #     self.mode.eval("save;")
    #     print('Sweep complete')
    #     neff = np.real(self.mode.getv("neff"))

    #     is_TE = abs(self.mode.getv("is_TE"))

    #     if self.close:
    #         self.close_instance()

    #     return neff, is_TE
