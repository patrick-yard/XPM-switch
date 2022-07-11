import numpy as np
from mode_solutions import Solver
import matplotlib.pyplot as plt
import os,json

class Waveguide(Solver):

    '''Class to simulate waveguide geometry for XPM.
        Parses input parameters for Solver and has functions to get and plot desired properties
        '''

    def __init__(self,material_parameters:dict,simulation_parameters:dict,mode:str,hide:bool,close:bool):
        
        '''Inputs:
            material_parameters: dictionary containing the waveguide dimensions and core/cladding materials
            simulation_parameters: dictionary containing paramters relating to the simulation region
            mode: string indicating which transverse mode type should be returned TE or TM
            hide: Bool deciding whether the GUI is hidden
            close: Bool deciding whether the GUI is closed once the simulation is finished'''

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
    

class DummyWaveguide():

    '''Dummy waveguide class, loads presaved data'''

    def __init__(self,folder):
        
        '''input:
        Folder should contain only one simulation data JSON file at a time.'''

        self.folder = folder

        self.load_data()
        self.combine_E_field
    
    def load_data(self) -> None:

        '''Load JSON file into dictionary'''

        self.file = [os.path.join(self.folder,f) for f in os.listdir(self.folder) if f.endswith(".JSON")][0]

        with open(self.file,'r') as jsonfile:
            self.data_dict = json.loads(jsonfile)

        return

    def combine_E_field(self) -> None:
        
        '''E field profiles split into real and imaginary parts (JSON doesnt like complex numbers apparently). This function recombines'''
        
        E_field_temp = []
        
        for E_real,E_imag in zip(self.data_dict['E_field_real'],self.data_dict['E_field_imag']):

            E_field_temp.append(np.array(E_real) + 1j*np.array(E_imag))
        
        self.data_dict['E_field'] = E_field_temp

        return
