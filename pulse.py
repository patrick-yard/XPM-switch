import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c,epsilon_0
from scipy.interpolate import interp2d

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

            self.waveguide.get_mode_profile(plot = False)
        
        self.mode_profile = self.waveguide.E_field[self.mode_idx]
        self.interp_real = interp2d(self.waveguide.y_array,self.waveguide.x_array, np.real( self.mode_profile))
        self.interp_imag = interp2d(self.waveguide.y_array,self.waveguide.x_array, np.imag( self.mode_profile))
        self.get_Aeff()
    
    def interpolated_E_field(self,x:int,y:int) -> complex:

        '''return interpolated E field at given point in FDE region '''

        return self.interp_real(y,x) + 1j*self.interp_imag(y,x)

    def get_interpolated_fields(self,num_points= 1001) -> tuple:

        '''Get interpolated E field across waveguide and whole FDE region
        inputs:
            num_points: number of interpolated points across entire FDE region
        returns
            interp_list: tuple containing interpolated E fields and x and y axes'''

        x_res = (max(self.waveguide.x_array)-min(self.waveguide.x_array))/num_points #defines resolution of interpolated field
        y_res = (max(self.waveguide.y_array)-min(self.waveguide.y_array))/num_points

        x_int = np.linspace(min(self.waveguide.x_array),max(self.waveguide.x_array),num_points) # x and y arrays for whole interpolated FDE area
        y_int = np.linspace(min(self.waveguide.y_array),max(self.waveguide.y_array),num_points)

        x_int_wg = np.arange(-self.waveguide.material_parameters['width']/2,self.waveguide.material_parameters['width']/2,x_res)# x and y arrays for just the waveguide
        y_int_wg = np.arange(-self.waveguide.material_parameters['height']/2,self.waveguide.material_parameters['height']/2,y_res)

        E_field_wg = self.interpolated_E_field(x_int_wg,y_int_wg)
        E_field = self.interpolated_E_field(x_int,y_int)

        self.interp_list = ((x_int,y_int),
                       (x_int_wg,y_int_wg),
                       E_field,
                       E_field_wg
                       )
        
        return self.interp_list

    
    def get_Aeff(self):

        int_axs,_,E_field,_ = self.get_interpolated_fields()
        
        denomenator = np.trapz(np.trapz(abs(E_field)**2 * abs(E_field)**2,x=int_axs[1]),x=int_axs[0])
        
        numerator = np.trapz(np.trapz(abs(E_field)**2,x=int_axs[1]),x=int_axs[0])**2

        self.Aeff = numerator/denomenator
        

         
    def get_vg(self):
        
        '''Function to find group velocity of given mode at lambda0'''
        
        if not hasattr(self.waveguide, 'vg'):

            self.waveguide.get_dispersion_info(plot = False)
        
        
        lambd0_idx = np.argmin([abs(wav - self.lambda0) for wav in self.waveguide.wav])
        self.f = self.waveguide.f[lambd0_idx]
        self.omega = 2*np.pi*self.f
        self.vg = self.waveguide.vg[self.mode_idx][lambd0_idx]

       

    def get_neff(self):
        
        '''Function to find effective index of given mode at lambda0'''
        
        if not hasattr(self.waveguide, 'neff'):

            self.waveguide.get_dispersion_info(plot = False)
        
        
        lambd0_idx = np.argmin([abs(wav - self.lambda0) for wav in self.waveguide.wav])

        self.neff = np.real(self.waveguide.neff[self.mode_idx][lambd0_idx])
        

    def get_peak_power(self):
        '''get peak power given pulse parameters'''

        self.peak_power = (self.power*1e-3)/(self.dt*1e-12) 

    def E(self,z:float,t:float) -> float:
            
        ''''Pulse E field as a function of distance along waveguide, z, and time ,t.'''

        A0 = 1#np.sqrt(4/(c*self.n*epsilon_0*self.Aeff)*self.peak_power)

        A = A0*np.exp(-((t-z/self.vg)/self.dt)**2)*np.exp(-((z-t*self.vg)/(self.dt*self.vg))**2) #tidy up when working

        return A*np.exp(1j*((2*np.pi/self.lambda0)*z-self.omega*t))

    def I(self,z:float,t:float) -> float:

        ''''Pulse intensity as a function of distance along waveguide, z, and time ,t.'''
        # (c*epsilon_0*self.n/2)
        return 1 * self.E(z,t) * self.E(z,t).conj()

    

        
            

