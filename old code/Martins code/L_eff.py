#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 22:24:14 2021

@author: Bielak Martin - bielak@optics.upol.cz
"""
import scipy.integrate as integrate
import scipy.constants as const

import numpy as np

###################
#   definitions   #
###################

k = lambda lamb: 2*np.pi /lamb
omega = lambda lamb: 2*np.pi * f(lamb)
f = lambda lamb: c/lamb
z_0 = lambda tau_p,v_g: v_g*tau_p
betha = lambda v_g: 1/v_g
#lambda_n = lambda lamb,n: lamb/n
#v_p = lambda n: c/n

# import a few constants
c = const.c
eps0 = const.epsilon_0

# pulse def.
A = lambda z,t,tau_p,v_g: np.exp(-((t-z/v_g)/tau_p)**2)*np.exp(-((z-t*v_g)/z_0(tau_p,v_g))**2)

def E(z,t,A0,tau_p,v_g,lamb):
    E_ = A0*A(z,t,tau_p,v_g)*np.exp(1j*(k(lamb)*z-omega(lamb)*t))
    return (E_)

def I(z,t,A0,tau_p,v_g,lamb,n):
    I_ = (c*eps0*n/2)*(E(z,t,A0,tau_p,v_g,lamb)*E(z,t,A0,tau_p,v_g,lamb).conjugate()).real
    return (I_)

# calculate power
P_W = lambda P_dBm: 1/1000 *10**(P_dBm/10) # convert dBm to W
P_peak = lambda P_avg, tau_p, f_rep:0.94*P_avg/(f_rep*tau_p)
A0 = lambda Aeff,n,P_avg,tau_p,f_rep: np.sqrt(4/(c*n*eps0*Aeff)*P_peak(P_avg,tau_p,f_rep))

#Gamma
dBetha = lambda v_g1, v_g2: betha(v_g1)-betha(v_g2)
gamma = lambda v_g1, v_g2, tau_p01, tau_p02, L: (0.44*2*np.pi)/(2*np.sqrt(2*np.log(2)))*1/np.sqrt(tau_p01**2+tau_p02**2)*np.abs(dBetha(v_g1,v_g2))*L
#for gamma >> 1
Leff_gamma = lambda v_g1, v_g2, tau_p01, tau_p02, L: np.sqrt(2)/gamma(v_g1, v_g2, tau_p01, tau_p02, L)*L

def L_eff(L,f_rep,Aeff01,tau_p01,v_g01,lamb01,Aeff02,tau_p02,v_g02,lamb02,P_avg01 = -27,P_avg02 = 10,n01 = 1.984675,n02 = 1.984675):
    L = L*10**(-3)
    f_rep = f_rep *10**6
    tau_p01 = tau_p01*10**(-12)
    tau_p02 = tau_p02*10**(-12)
    tau_delay = 1/2*((L/v_g01)-(L/v_g02))
    Aeff01 = Aeff01 *10**(-12)
    Aeff02 = Aeff02 *10**(-12)
    A01 = A0(Aeff01,n01,P_W(P_avg01),tau_p01,f_rep)
    A02 = A0(Aeff02,n02,P_W(P_avg02),tau_p02,f_rep)
    lamb01 = lamb01*10**(-9)
    lamb02 = lamb02*10**(-9)
    
    #integration bounds (substitution of infinity in dt)
    endT = 2*L/v_g01 #L/v_g01 = time of propagation of the 1. pulse through the whole waveguide
    tau_p_ = min(tau_p01,tau_p02)
    
    # denominator
    t_array = np.linspace(0,endT,3000)
    I_temp = lambda t: I(L/2,t,A01,tau_p01,v_g01,lamb01,n01)*I(L/2,t-tau_delay,A02,tau_p02,v_g02,lamb02,n02)
    dt = t_array[1]
    int_array = I_temp(t_array)*dt
    denom = np.sum(int_array)
    
    # numerator
    t_array = np.linspace(0,endT,5000)
    z_array = np.linspace(0,L,5000)
    I_temp = lambda t,z: I(z,t,A01,tau_p01,v_g01,lamb01,n01)*I(z,t-tau_delay,A02,tau_p02,v_g02,lamb02,n02)
    res= np.zeros((len(t_array),len(z_array)))
    for i,t in enumerate(t_array):
        res[i] = I_temp(t,z_array)
    int_res = np.trapz(np.trapz(res,x = t_array),x = z_array)
    
    l_eff = int_res/denom
    return (l_eff)

def Print_info(L,f_rep,Aeff01,tau_p01,v_g01,lamb01,Aeff02,tau_p02,v_g02,lamb02,P_avg01 = -27,P_avg02 = 10,n01 = 1.984675,n02 = 1.984675):
    l_eff = L_eff(L,f_rep,Aeff01,tau_p01,v_g01,lamb01,Aeff02,tau_p02,v_g02,lamb02,P_avg01,P_avg02,n01, n02)
    
    L = L*10**(-3)
    f_rep = f_rep *10**6
    tau_p01 = tau_p01*10**(-12)
    tau_p02 = tau_p02*10**(-12)
    Aeff01 = Aeff01 *10**(-12)
    Aeff02 = Aeff02 *10**(-12)
    lamb01 = lamb01*10**(-9)
    lamb02 = lamb02*10**(-9)
    
    gamma_vall = gamma(v_g01, v_g02, tau_p01, tau_p02, L)
    leff_gamma = Leff_gamma(v_g01, v_g02, tau_p01, tau_p02, L)
    print ('Gamma = \t\t' + str(gamma_vall))
    print ('\nFor Gamma >> 1:')
    print ('Leff_gamma = \t\t' + str(leff_gamma))
    print ('\nFor Gamma << 1:')
    print ('Leff = L = \t\t' + str(L))
    print ('\n"True" calculation')
    print ('Len. of waveguide = \t' + str(L) + ' m')
    print ('Interaction len = \t' + str(l_eff) + ' m')
    return ([gamma_vall,l_eff,leff_gamma,L])


if __name__ == '__main__':
    L = 30 #length of waveguide in mm
    f_rep = 50 #MHz
    
    #Pulse 1 - sig.: TE0, 1555nm, 6.18ps
    Aeff01 = 1.258795448537488 #um^2
    tau_p01 = 6.18 #pulse length in ps
    v_g01 = 152228482.56849426 #propagation speed in [m/s]
    lamb01 = 1555 #lamb in nm
    
    #Pulse 2 - pump: TE1, 1560nm, 3.87ps
    Aeff02 = 1.6575927059408382 #um^2
    tau_p02 = 3.87 #pulse length in ps
    v_g02 = 148663622.78341362 #propagation speed in [m/s]
    lamb02 = 1560 #lamb in nm
    
    l_eff = L_eff(L,f_rep,Aeff01,tau_p01,v_g01,lamb01,Aeff02,tau_p02,v_g02,lamb02)
    print (l_eff)
    Print_info(L,f_rep,Aeff01,tau_p01,v_g01,lamb01,Aeff02,tau_p02,v_g02,lamb02)