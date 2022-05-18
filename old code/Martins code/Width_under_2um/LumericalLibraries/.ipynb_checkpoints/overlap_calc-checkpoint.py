import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
from LumericalLibraries.mode_solutions import mode_sol

class Overlap_calc:

	'''Class to calculate various overlap integrals '''

	def __init__(self,sim_parameters,modes,initialise = True):


		N_points_Disp = 21
		start_wavelength = 1.555e-6
		stop_wavelength = 1.560e-6
		d_wavelength = (stop_wavelength-start_wavelength) / (N_points_Disp-1)
		wavelength = 1.55e-6
		Lc = 1e-6
		substrate_height = 2e-6
		num_TE = 3
		mesh_dx = 0.05e-6
		mesh_dy = 0.05e-6


		self.simulation_parameters = {
		"wavelength" : wavelength,
		"start_wavelength":start_wavelength,
		"stop_wavelength":stop_wavelength,
		"N_points_Disp":N_points_Disp,
		"d_wavelength" : d_wavelength,
		"Lc" : Lc,
		"N_points_Disp" : N_points_Disp,
		"mesh_dx" : mesh_dx,
		"mesh_dy" : mesh_dy,
		"substrate_height":substrate_height,
		"num_TE" : num_TE
		}

		
		self.simulation_parameters.update(sim_parameters) #paramaters to be passed to lumerical
		if initialise:
			self.MODE = mode_sol(self.simulation_parameters)
		self.mode_idx = modes

		self.x_array = 0
		self.y_array = 0
		self.E_fields = 0
		

	def Get_E_field(self,freq_sweep = False,plot = False):

		""" Extract E fields from mode solutions option to plot fields and refractive index profile"""

		self.E_fields,self.index_profile,self.x_array,self.y_array,self.vgs = self.MODE.Extract_E_field(freq_sweep=freq_sweep)
		
		
		if plot:
			self.plot_fields()

		return self.E_fields,self.index_profile,self.x_array,self.y_array,self.vgs

	def interpolate_field(self):
		
		'''Interpolate E field -> use if mesh size is not uniform across FDE region'''
		
		int_E_field = []
	
		for E in  self.E_fields:
			E_real = np.real(E)
			E_imag = np.imag(E)
			int_field_real = interp2d(self.y_array,self.x_array, E_real)
			int_field_imag = interp2d(self.y_array,self.x_array, E_imag)
			int_E_field.append([int_field_real,int_field_imag])
		return int_E_field

	def interpolate_index(self):
		
		'''Interpolate index profile -> use if mesh size is not uniform across FDE region'''
		

		return interp2d(self.x_array,self.y_array, self.index_profile)

	def plot_index(self,profile):
		
		fig = plt.figure()
		fig, ax = plt.subplots()
		plot = ax.contourf(self.x_array*1e6,self.y_array*1e6,index_profile,levels = 500)
		cbar = fig.colorbar(plot,label = 'index_profile')
		
		ax.set_xlabel('width (x) (um)')
		ax.set_ylabel('height (y) (um)')
		plt.show()

		return
	def plot_fields(self,index_profile = False):
        
		for E in  self.E_fields:


			fig = plt.figure()
			fig, ax1 = plt.subplots()
			plot = ax1.contourf(self.x_array*1e6,self.y_array*1e6,np.real(E.T),levels = 500)
			cbar = fig.colorbar(plot,label = 'Ex')
			
			ax1.set_title('Fields from Lumerical')
			ax1.set_xlabel('width (x) (um)')
			ax1.set_ylabel('height (y) (um)')
			plt.show()
		

		if index_profile:
			self.plot_index()


		return
	

	def plot_int_fields(self,interpolated_Efield,int_index_profile,x,y):


		plt.rcParams.update({'font.size':18})


		for i,i_E in enumerate(interpolated_Efield):

			E_real = i_E[0](y,x).T
		
			# fig = plt.figure()
			fig, ax = plt.subplots(figsize = (6,3))
			plot = ax.contourf(x*1e6,y*1e6,abs(E_real)**2,levels = 500)
			cbar = fig.colorbar(plot,ticks = [0,0.5,1],label = r'$\vert E_x \vert^2$')
			
			ax.set_xlabel(r'width (x) ($\mu$m)')
			ax.set_ylabel(r'height (y) ($\mu$m)')
			# ax.set_title('Interpolated fields')
			# fig.savefig(f'mode_profile_{i}.png',bbox_inches = 'tight')
			plt.show()
		fig = plt.figure()
		fig, ax = plt.subplots()
		int_profile = int_index_profile(x,y)
		
		plot = ax.contourf(x*1e6,y*1e6,int_profile,levels = 500)
		cbar = fig.colorbar(plot,label = 'index_profile')
		
		ax.set_xlabel('width (x) (um)')
		ax.set_ylabel('height (y) (um)')
		plt.show()
		return

	def interp_axes(self,num_int_points):

		x_res_int = (max(self.x_array)-min(self.x_array))/num_int_points #defines resolution of interpolated field
		y_res_int = (max(self.y_array)-min(self.y_array))/num_int_points

		x_int = np.linspace(min(self.x_array),max(self.x_array),num_int_points) # x and y arrays for whole interpolated FDE area
		y_int = np.linspace(min(self.y_array),max(self.y_array),num_int_points)
	    
		x_int_wg = np.arange(-self.simulation_parameters['width']/2,self.simulation_parameters['width']/2,x_res_int)# x and y arrays for just the waveguide
		y_int_wg = np.arange(-self.simulation_parameters['height']/2,self.simulation_parameters['height']/2,y_res_int)

		return x_int,y_int,x_int_wg,y_int_wg

	def overlap_calc(self,modes,num_int_points,int_E_field):

		'''Nonlinar overlap : doesnt include refractive index in denomentaor -> Majority of mode is in core TE'''

		x_int,y_int,x_int_wg,y_int_wg = self.interp_axes(num_int_points)

		field1 = int_E_field[modes[0]]
		field2 = int_E_field[modes[1]]

		E_field_wg1 = field1[0](y_int_wg,x_int_wg) + 1j*field1[1](y_int_wg,x_int_wg)
		E_field_wg2 = field2[0](y_int_wg,x_int_wg) + 1j*field2[1](y_int_wg,x_int_wg)
		
		E_field_tot1 = field1[0](y_int,x_int) + 1j*field1[1](y_int,x_int)
		E_field_tot2 = field2[0](y_int,x_int) + 1j*field2[1](y_int,x_int)

		numerator = np.trapz(np.trapz(abs(E_field_wg1)**2 * abs(E_field_wg2)**2,x=y_int_wg),x=x_int_wg)
		denomenator = np.trapz(np.trapz(abs(E_field_tot1)**2,x=y_int),x=x_int)*np.trapz(np.trapz(abs(E_field_tot2)**2,x=y_int),x=x_int)

		overlap = numerator/denomenator

		return overlap 

	def linear_overlap(self,modes,num_int_points,int_E_field):

		'''Linear overlap -- mainly for testing'''

		x_int, y_int, _ , _ = self.interp_axes(num_int_points)

		field1 = int_E_field[modes[0]]
		field2 = int_E_field[modes[1]]

		
		E_field_tot1 = field1[0](y_int,x_int) + 1j*field1[1](y_int,x_int)
		E_field_tot2 = field2[0](y_int,x_int) + 1j*field2[1](y_int,x_int)

		numerator = np.trapz(np.trapz(E_field_tot1.conj()*E_field_tot2,x=y_int),x=x_int)
		denomenator = np.trapz(np.trapz(abs(E_field_tot1)**2,x=y_int),x=x_int)*np.trapz(np.trapz(abs(E_field_tot2)**2,x=y_int),x=x_int)

		overlap = abs(numerator)**2/denomenator

		return overlap 
 
	def NL_overlap(self,freq_sweep = 0,num_int_points = 1001):

		'''Calculates the nonlinear overlap between different transverse field modes'''

		self.Get_E_field(freq_sweep=freq_sweep)
		
		int_E_field = np.array(self.interpolate_field())
		
		int_n = self.interpolate_index()
		
		x_int, y_int, _ , _ = self.interp_axes(num_int_points)
		
		nonlinear_overlap = self.overlap_calc(self.mode_idx,num_int_points,int_E_field) ### desired overlap
		
		lin_ols = [self.linear_overlap([0,0],num_int_points,int_E_field),self.linear_overlap([1,1],num_int_points,int_E_field),self.linear_overlap([0,1],num_int_points,int_E_field)] ### desired overlap
		
		self.plot_int_fields(int_E_field,int_n,x_int,y_int)
		
		TE0_overlap = self.overlap_calc([0,0],num_int_points,int_E_field) ### overlap for TE0 with itself for testing
		
		TE1_overlap = self.overlap_calc([1,1],num_int_points,int_E_field)
		

		return nonlinear_overlap,TE0_overlap,TE1_overlap,lin_ols
	
