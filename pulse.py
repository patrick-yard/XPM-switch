import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c,epsilon_0


class Pulse():

    '''Pulse object. Inputs:
    mode_idx: which mode should be associated to this pulse? TE/TM is set in waveguide class
    lambd0: pulse central wavelength in nm
    dt: pulsewidth in ps
    power: power per pulse in mW - average_power/rep_rate
    waveguide: Waveguide class containing geometry and lumerical simulation components
    '''

    def __init__(self,
                mode_idx:int,
                lambda0:float,
                dt:float,
                power:float,
                waveguide:object
                ):

        self.mode_idx = mode_idx
        self.lambda0 = lambda0
        self.dt = dt
        self.power = power
        self.peak_power = self.get_peak_power()
        self.waveguide = waveguide
        self.n = self.waveguide.n
        self.get_mode_profile()
        self.get_vg()
        self.get_neff()
        self.get_peak_power()
    
        
    def get_mode_profile(self):
        
        '''Function to pick out field profile of the mode specificed by self.mode_idx'''

        if not hasattr(self.waveguide, 'E_field'):

            self.get_mode_profile(plot = False)
        
        self.mode_profile = self.waveguide.E_field[self.mode_idx]
        self.aeff = [] # need to integrate mode_profile
    
    def get_vg(self):
        
        '''Function to find group velocity of given mode at lambda0'''
        
        if not hasattr(self.waveguide, 'vg'):

            self.get_dispersion_info(plot = False)
        
        
        lambd0_idx = np.argmin([abs(wav - self.lambda0) for wav in self.waveguide.wav])
        self.f = self.waveguide.f[lambd0_idx]
        self.omega = 2*np.pi*self.f
        self.vg = self.waveguide.vg[self.mode_idx][lambd0_idx]


    def get_neff(self):
        
        '''Function to find effective index of given mode at lambda0'''
        
        if not hasattr(self.waveguide, 'neff'):

            self.get_dispersion_info(plot = False)
        
        
        lambd0_idx = np.argmin([abs(wav - self.lambda0) for wav in self.waveguide.wav])

        self.neff = self.waveguide.neff[self.mode_idx][lambd0_idx]

    def get_peak_power(self):
        '''get peak power given pulse parameters'''

        self.peak_power = (self.power*1e-3)/(self.dt*1e-12) 

    def E(self,z:float,t:float) -> float:
            
        ''''Pulse E field as a function of distance along waveguide, z, and time ,t.'''

        A0 = np.sqrt(4/(c*self.n*epsilon_0*self.Aeff)*self.peak_power)

        A = A0*np.exp(-((t-z/self.v_g)/self.dt)**2)*np.exp(-((z-t*self.v_g)/(self.dt*self.v_g))**2) #tidy up when working

        return A*np.exp(1j*(self.neff*z-self.omega*t))

    def I(self,z:float,t:float) -> float:

        ''''Pulse intensity as a function of distance along waveguide, z, and time ,t.'''

        return (c*epsilon_0*self.n/2) * self.E(z,t) * self.E(z,t).conj()

    

        
            

