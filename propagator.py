from lib2to3.pgen2.token import LEFTSHIFTEQUAL
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

        T_slow = 2*self.length/min([self.pump.vg,self.signal.vg]) #time for slowest pulse to propagate whole waveguide

        t_array = np.linspace(0,5*T_slow,5000)
        z_array = np.linspace(0,self.length,5000)

        num_integrand = self.pump.I(z_array,t_array)*self.signal.I(z_array,t_array)
        denom_integrand = self.pump.I(self.length/2,t_array)*self.signal.I(self.length/2,t_array)

        numerator = np.trapz(np.trapz(num_integrand,x = t_array),x = z_array)
        denomenator = np.trapz(denom_integrand,x = t_array)
        
        return numerator/denomenator
    
    def get_wg_indices(self) -> tuple:

        '''Find field profile indices corresponding to the waveguide location '''

        width_start = -self.waveguide.material_parameters['width']/2
        width_start_idx = np.argmin([w - width_start for w in self.waveguide.x_array])

        width_end = self.waveguide.material_parameters['width']/2
        width_end_idx = np.argmin([w - width_end for w in self.waveguide.x_array])

        height_start = -self.waveguide.material_parameters['height']/2
        height_start_idx = np.argmin([w - height_start for w in self.waveguide.x_array])

        height_end = self.waveguide.material_parameters['height']/2
        height_end_idx = np.argmin([w - height_end for w in self.waveguide.y_array])

        return (width_start_idx,width_end_idx,height_start_idx,height_end_idx)



        return
    def calculate_overlap(self,linear = False) -> float:

        '''Calcualte nonlinear overlap between pulse profiles
            Input:
                linear: If True calculates linear (field) overlap, if False calculates nonlinear (intesity) overlap
        '''

        pump_field = self.pump.mode_profile
        signal_field = self.signal.mode_profile

        width_start,width_end,height_start,height_end = self.get_wg_indices()

        pump_field_wg = pump_field[width_start:width_end,height_start:height_end]
        signal_field_wg = signal_field[width_start:width_end,height_start:height_end]

        wg_height_array = self.waveguide.yarray[height_start:height_end]
        wg_width_array = self.waveguide.xarray[width_start:width_end]
        
        if linear:
            numerator = np.trapz(pump_field.conj().T*signal_field,x=wg_height_array,x=wg_width_array)
        else:
            numerator = np.trapz(np.trapz(abs(pump_field_wg)**2 * abs(signal_field_wg)**2,x=wg_height_array),x=wg_width_array)
        
        denomenator = np.trapz(np.trapz(abs(pump_field_wg)**2,x=self.waveguide.yarray),x=self.waveguide.xarray)*np.trapz(np.trapz(abs(signal_field_wg)**2,x=self.waveguide.yarray),x=self.waveguide.xarray)

        return numerator/denomenator

    def calculate_phase(self) -> float:

        '''Calculate XPM phase shift imparted by pump onto signal'''

        prefactor = (4*np.pi*self.n2)/self.signal.wav
        
        ref_indices = (self.pump.ng*self.signal.ng)/self.Waveguide.n

        overlap = self.calculate_overlap(linear = False)

        power = self.pump.peak_power

        Leff = self.caculate_Leff()
        
        phase = prefactor * ref_indices *  overlap * power * Leff

        return phase