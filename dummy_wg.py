
import numpy as np
import os,json
 
class DummyWaveguide():

    '''Dummy waveguide class, loads presaved data'''

    def __init__(self,folder):
        
        '''input:
        Folder should contain only one simulation data JSON file at a time.'''

        self.folder = folder

        self.load_data()
        self.combine_E_field()
        # print(self.data_dict)
        self.__dict__ = self.data_dict
    
    def load_data(self) -> None:

        '''Load JSON file into dictionary'''
        
        self.file = [os.path.join(self.folder,f) for f in os.listdir(self.folder) if f.endswith(".JSON")][0]

        with open(self.file,'r') as jsonfile:
            self.data_dict = json.loads(json.load(jsonfile))

        return

    def combine_E_field(self) -> None:
        
        '''E field profiles split into real and imaginary parts (JSON doesnt like complex numbers apparently). This function recombines'''
        
        E_field_temp = []
        
        for E_real,E_imag in zip(self.data_dict['E_field_real'],self.data_dict['E_field_imag']):

            E_field_temp.append(np.array(E_real) + 1j*np.array(E_imag))
        
        self.data_dict['E_field'] = E_field_temp

        return
