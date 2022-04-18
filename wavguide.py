import numpy as np
import os
import sys
from scipy.constants import c
from mode_solutions import Solver

GENERAL_SIM_PARAMS = {}

class Waveguide(Solver):

    def __init__(self,
                waveguide_parameters:dict,
                general_params = GENERAL_SIM_PARAMS,
                mode = 'TE',
                hide = True,
                close = True
                ):

        '''Class representing a Waveguide component. Inputs waveguide parameters and simulation parameters'''

        self.waveguide_parameters = waveguide_parameters
        self.waveguide_parameters.update(general_params)

        super.__init__(self.waveguide_parameters,mode = mode, hide = hide, close = close)
    
    def get_field_profile(self,index:int or list) -> dict():

        '''Function to get field profile for specificed mode index or indices'''

        return

    def get_dispersion_information(self) -> dict():

        '''Function to return simulated n_eff vs f'''

        return

    def parameter_sweep(self,param_dict:dict)-> dict:

        '''Function to sweep a given waveguide parameter'''

        return