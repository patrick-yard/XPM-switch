from itertools import product
from lib2to3.pgen2.token import LEFTSHIFTEQUAL
from signal import signal
import numpy as np
import matplotlib.pyplot as plt
from pulse import Pulse


class Propagator():

    '''Propagtor object with functionality to co-propagate two pulses and determine phase shifts
    Inputs:
        Waveguide: Wavguide class defines waveguide properties and dispersion info
        pump_params: dict containing pulse parameters for pump pulse, needs to be correct format to be passed directly to Pulse class
        signal_params: dict containing pulse parameters for signal pulse
        length: length to propagate over in mm
        n0: Material refractive index
        n2: nonlinear refractive index

    '''

    def __init__(self,
                Waveguide:object,
                pump_params:dict,
                signal_params:dict,
                length:float,
                n2:float,
                ):

        self.Waveguide = Waveguide
        self.pump_params = pump_params
        self.signal_params = signal_params
        self.length = length
        self.n2 = n2

        self.pump = Pulse(waveguide = self.Waveguide,**self.pump_params)
        self.signal = Pulse(waveguide = self.Waveguide,**self.signal_params)
    
    def caculate_Leff(self) -> float:

        '''Calculates the effective length of two propagating pulses'''

        T_slow = self.length/min([self.pump.vg,self.signal.vg]) #time for slowest pulse to propagate whole waveguide
        tau = 0.5*self.length*(1/self.pump.vg - 1/self.signal.vg)
        
        t_array = np.linspace(0,2*T_slow,2000)
        z_array = np.linspace(0,self.length,2000)

        num_integrand = np.array([self.pump.I(z_array,t)*self.signal.I(z_array,t-tau) for t in t_array])
        denom_integrand =self.pump.I(self.length/2,t_array)*self.signal.I(self.length/2,t_array-tau)

        numerator = np.trapz(np.trapz(num_integrand,x = t_array),x = z_array)
        denomenator = np.trapz(denom_integrand,x = t_array)
        
        return numerator/denomenator
    
    def get_maximum_leff(self) -> float:
        
        '''Returns effective length in the case of full walkoff'''

        prefactor  = (4*np.sqrt(np.log(2)))/(0.441271*np.pi*2)

        numerator = prefactor * np.sqrt(self.pump.dt**2 + self.signal.dt**2)

        denomenator = abs((1/self.pump.vg - 1/self.signal.vg))

        return numerator/denomenator

    def calculate_overlap(self,linear: bool = False) -> float:

        '''Calcualte nonlinear overlap between pulse profiles
            Input:
                linear: If True calculates linear (field) overlap, if False calculates nonlinear (intesity) overlap
        '''

        pump_field = self.pump.interp_list[2]
        signal_field = self.signal.interp_list[2]
        
        x,y = self.pump.interp_list[0]
        
        pump_field_wg = self.pump.interp_list[3]
        signal_field_wg = self.pump.interp_list[3]

        x_wg,y_wg = self.pump.interp_list[1]
        print(len(x_wg),len(y_wg))
        print(pump_field_wg.shape,signal_field_wg.shape)
        if linear:
                       
            numerator = abs(np.trapz(np.trapz(pump_field.conj()*signal_field,x=y),x=x))**2
        else:
            # numerator = np.trapz(np.trapz(abs(pump_field_wg)**2 * abs(signal_field_wg)**2,x=y_wg),x=x_wg)
            numerator = np.trapz(np.trapz(abs(pump_field)**2 * abs(signal_field)**2,x=y),x=x)
        
        denomenator = np.trapz(np.trapz(abs(pump_field)**2,x=y),x=x)*np.trapz(np.trapz(abs(signal_field)**2,x=y),x=x)
        # print(numerator,denomenator)
        return numerator/denomenator

    def calculate_phase(self) -> float:

        '''Calculate XPM phase shift imparted by pump onto signal'''

        prefactor = (4*np.pi*self.n2)/self.signal.wav
        
        ref_indices = (self.pump.ng*self.signal.ng)/self.Waveguide.n**2

        overlap = self.calculate_overlap(linear = False)

        power = self.pump.peak_power

        Leff = self.caculate_Leff()

        phase = prefactor * ref_indices *  overlap * power * Leff

        return phase




          # def get_wg_indices(self) -> tuple:

    #     '''Find field profile indices corresponding to the waveguide location '''

    #     width_start = -self.Waveguide.material_parameters['width']/2
        
    #     width_start_idx = np.argmin([abs(w - width_start) for w in self.Waveguide.x_array])
       
    #     width_end = self.Waveguide.material_parameters['width']/2
    #     width_end_idx = np.argmin([abs(w - width_end) for w in self.Waveguide.x_array])

    #     height_start = -self.Waveguide.material_parameters['height']/2
    #     height_start_idx = np.argmin([abs(w - height_start) for w in self.Waveguide.y_array])

    #     height_end = self.Waveguide.material_parameters['height']/2
    #     height_end_idx = np.argmin([abs(w - height_end) for w in self.Waveguide.y_array])

    #     return (width_start_idx,width_end_idx,height_start_idx,height_end_idx)




    # def interpolate_fields(self) -> tuple:

    #     '''Function to interpolate electric field profile and corresponding axes'''


    #     return
    # def calculate_overlap(self,linear: bool = False) -> float:

    #     '''Calcualte nonlinear overlap between pulse profiles
    #         Input:
    #             linear: If True calculates linear (field) overlap, if False calculates nonlinear (intesity) overlap
    #     '''

    #     pump_field = self.pump.mode_profile
    #     signal_field = self.signal.mode_profile


    #     width_start,width_end,height_start,height_end = self.get_wg_indices()
    #     # print(width_start,width_end,height_start,height_end)
        
        
    #     pump_field_wg = pump_field[width_start:width_end,height_start:height_end]
    #     signal_field_wg = signal_field[width_start:width_end,height_start:height_end]

    #     wg_height_array = self.Waveguide.y_array[height_start:height_end]
    #     wg_width_array = self.Waveguide.x_array[width_start:width_end]


    #     # print(self.Waveguide.y_array[height_start],self.Waveguide.y_array[height_end])
    #     # print(self.Waveguide.x_array[width_start],self.Waveguide.x_array[width_end])
    #     # print(wg_height_array.shape,wg_width_array.shape)
    #     # print(pump_field_wg.shape,signal_field_wg.shape)
               
       
    #     if linear:
                       
    #         numerator = abs(np.trapz(np.trapz(pump_field.conj()*signal_field,x=self.Waveguide.y_array),x=self.Waveguide.x_array))**2
    #     else:
    #         numerator = np.trapz(np.trapz(abs(pump_field_wg)**2 * abs(signal_field_wg)**2,x=wg_height_array),x=wg_width_array)
    #         # numerator = np.trapz(np.trapz(abs(pump_field)**2 * abs(signal_field)**2,x=self.Waveguide.y_array),x=self.Waveguide.x_array)
        
    #     denomenator = np.trapz(np.trapz(abs(pump_field)**2,x=self.Waveguide.y_array),x=self.Waveguide.x_array)*np.trapz(np.trapz(abs(signal_field)**2,x=self.Waveguide.y_array),x=self.Waveguide.x_array)
    #     # print(numerator,denomenator)
    #     return numerator/denomenator