import numpy as np
from mode_solutions import Solver
import matplotlib.pyplot as plt

class Waveguide(Solver):

    '''Class to simulate waveguide geometry for XPM.
        Parses input parameters for Solver and has functions to get and plot desired properties'''

    def __init__(self,material_parameters:dict,simulation_parameters:dict,mode:str,hide:bool,close:bool):

        self.material_parameters = material_parameters
        self.simulation_parameters = simulation_parameters
        self.other_params = (mode,hide,close)
        

        self.total_params = dict()
        self.total_params.update(self.material_parameters)
        self.total_params.update(self.simulation_parameters)
        self.n = material_parameters['n'] # TODO: read this from lumerical

        super().__init__(self.total_params,*self.other_params)


    def get_mode_profile(self, plot = True) -> np.array:
        
        ''' Return mode profiles for all guided modes and unpacks key properties '''

        field_dict = self.Extract_E_field()

        self.E_field = field_dict['E_field']
        self.index_profile = field_dict['index_profile']
        self.x_array =  field_dict['x_array']
        self.y_array =  field_dict['y_array']

        if plot:
            self.plot_field_profiles()
        
        return field_dict


    def plot_field_profiles(self):
        
        if not hasattr(self,'E_field'):

            self.get_mode_profile(plot = False)

        fig,axs = plt.subplots(nrows = 1, ncols = 1 + len(self.E_field),figsize = (15,5))

        for i,ax in enumerate(axs):

            if i == 0:
                plot_array = np.real(self.index_profile)
            else:
                plot_array = abs(self.E_field[i-1].T)**2
            
            ax.imshow(plot_array)
        
        plt.show()

        return
    def get_dispersion_info(self,plot = True) -> dict:

        '''get dispersion info. Returns dictionary of dispersion data and unpacks properties '''

        self.disp_dict = self.freq_sweep()

        self.vg = self.disp_dict['vg']
        self.D = self.disp_dict['D']
        self.neff = self.disp_dict['neff']
        self.wav = self.disp_dict['wav']
        self.f = self.disp_dict['f']

        if plot:
            self.plot_dispersion_data()

        return self.disp_dict

    def plot_dispersion_data(self):
        
        ''''Plots dispersion data'''

        if not hasattr(self,'disp_dict'):
            
             self.get_dispersion_info(plot = False)

        fig,axs = plt.subplots(nrows = 1, ncols = 3,figsize = (15,5))
        keys = ('neff','vg','D')
        for ax,key in zip(axs,keys):

            ax.set_ylabel(key)
            ax.set_xlabel('Wavelength (nm)')

            [ax.plot(self.wav[0]*1e9,data) for data in self.disp_dict[key]]

        fig.tight_layout()
        plt.show()

        return