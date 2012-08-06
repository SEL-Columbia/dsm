## Mitchell Lee
## Subfunctions for smartAllocate.py
## Script Begun August 5, 2012

from numpy import pi, sin, cos, tan, arccos, arcsin
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

# dates is an array 8760 x 6 columns 

# unvectorized

# this function takes the normal to sun, direct beam insolation from a
# solar database and converts it to insolation on a fixed panel
# accounting for both direct, ground reflected, and diffuse 

def resourceCalc(date, sigma, phi_c, I_B, lats, rho):
    '''
    inputs
    ------
    date - datetime object
    sigma - collector angle
    rho - factor for ground reflectance
    '''

    #time = datenum(dates);
    
    #sigma = (sigma)*pi/180  #collector angle sigma from ground
    sigma = sp.radians(sigma)
    L = lats
    L = sp.radians(L)
    phi_c = sp.radians(phi_c)

    #n = ceil((1:length(time))/24)';
    # give day of year
    n = int(date.strftime('%j'))


    delta = 23.45 * pi / 180. * sin(2 * pi / 365. * (n - 81))
    time_solar = int(date.strftime('%H'))+1
    H = 2 * pi / 24 * (12 - time_solar)
    
    beta = arcsin(cos(L) * cos(delta) * cos(H) + sin(L) * sin(delta))
    phi_s = arcsin(cos(delta) * sin(H) / cos(beta))
    
    # C ambient light correction
    C = 0.095 + 0.04 * sin(360 / 365. * (n - 100))
    # Stationary Collector

    I_C = I_B * (cos(beta) * cos(phi_s - phi_c) * sin(sigma)
                 + sin(beta) * cos(sigma)
                 + C * (1+cos(sigma)) / 2. + rho * (sin(beta)+C) * (1-cos(sigma))/2.)
    if I_C < 0:
        I_C = 0

    return I_C
	
# Basic Energy Balance for a Single Time Step
	
def energyBalance(batChar,supply,demand,batCap,batMin)  
	batCharNew = batChar + supply + np.sum(demand)
		supplyM = demand
		if batCharNew > batCap
			batCharNew = batCap
		elseif batCharNew < batMin
			batCharNew = batMin
			# Distributing remaining electricity
			# Electricity distributed proportionally to energy request shapes
			supplyM = supplyM/np.sum(supplyM)*(batChar-batMin)
	return batCharNew, supplyM
	
	