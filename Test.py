
%matplotlib inline
import numpy as np
import pylab as pl
import scipy.special as special
from scipy.integrate import quad

# Set plot parameters to make beautiful plots
pl.rcParams['figure.figsize']  = 12, 7.5
pl.rcParams['lines.linewidth'] = 1.5
pl.rcParams['font.family']     = 'serif'
pl.rcParams['font.weight']     = 'bold'
pl.rcParams['font.size']       = 20  
pl.rcParams['font.sans-serif'] = 'serif'
pl.rcParams['text.usetex']     = True
pl.rcParams['axes.linewidth']  = 1.5
pl.rcParams['axes.titlesize']  = 'medium'
pl.rcParams['axes.labelsize']  = 'large'

pl.rcParams['xtick.major.size'] = 8     
pl.rcParams['xtick.minor.size'] = 4     
pl.rcParams['xtick.major.pad']  = 8     
pl.rcParams['xtick.minor.pad']  = 8     
pl.rcParams['xtick.color']      = 'k'     
pl.rcParams['xtick.labelsize']  = 'large'
pl.rcParams['xtick.direction']  = 'in'    

pl.rcParams['ytick.major.size'] = 8     
pl.rcParams['ytick.minor.size'] = 4     
pl.rcParams['ytick.major.pad']  = 8     
pl.rcParams['ytick.minor.pad']  = 8     
pl.rcParams['ytick.color']      = 'k'     
pl.rcParams['ytick.labelsize']  = 'large'
pl.rcParams['ytick.direction']  = 'in'


# EVERYTHING IS 1 units

#constants
e        = 1.     #electron charge
m        = 1.     #electron mass
c        = 1.     #speed of light
epsilon0 = 1.     #permittivity of free space
epsilon  = -1.    #sign of electron charge

#parameters
B     = 1.         #background B strength
n_e   = 1.         #electron number density cm^-3
w_T   = 1.         #dimensionless electron temp. k_B T / m c^2
theta = np.pi / 3. #observer angle

#derived quantities
omega_p = np.sqrt(n_e * e**2. / (m * epsilon0))     # plasma frequency
omega_c = e * B / (m * c)                           # cyclotron frequency

#we first need the plasma dispersion function

def Z_integrand(xi, zeta):
    prefactor   = 1. / np.sqrt(np.pi)
    numerator   = np.exp(-xi**2.)
    
#    denominator = xi - zeta  #included in quad with weight type 'cauchy' passed to quad
    denominator = 1.
    
    return prefactor * numerator / denominator


#seems to work up to |zeta| = 625 where it's approx. -/+ 0.002 (so negative zeta yields +0.002)
def Z(zeta): 
    if(np.abs(zeta) < 2.):
        int_limit = 10. * np.abs(zeta)
    elif(np.abs(zeta) > 2. and np.abs(zeta) < 130.):
        int_limit = 2. * np.abs(zeta)
    else:
        int_limit = 1.5 * np.abs(zeta)
    
    imag_part = 1j * np.pi * Z_integrand(zeta, zeta)
    
    if(zeta != 0):
        ans = quad(lambda xi: Z_integrand(xi, zeta), -int_limit, int_limit, weight='cauchy', wvar=zeta)[0]
    else:
        ans = 0.
        
    return ans + imag_part

def K_12_summand(n, omega):
    k_perp = omega / c * np.sin(theta)                  # wavevector perp component n = 1 approximation
    k_z    = omega / c * np.cos(theta)                  # wavevector parallel comp. n = 1 approximation
    lambd   = k_perp**2. * w_T**2. / (2. * omega_c**2.) # just a defined parameter
    prefactor = 1j * omega_p**2. * np.exp(-lambd) / (omega * k_z * w_T)
    zeta = (omega + n * omega_c) / (k_z * w_T)
    
    if(np.abs(zeta) > 625):
        print 'zeta out of range of PDF'
    
    term1 = n * (special.iv(n, lambd) - special.ivp(n, lambd)) * Z(zeta)
    ans = prefactor * term1
    return ans

def K_12(terms, omega):
    ans = 0.
    for i in range(-terms, terms):
        ans += K_12_summand(i, omega)
#        print i, ans  
    return ans

#print K_12(number of sum terms evaluated, omega)

print K_12(100, 40.0)

#ACTUAL REASONABLE VALUES

#constants
e        = 4.80320680e-10    #electron charge
m        = 9.1093826e-28     #electron mass
c        = 2.99792458e10     #speed of light
epsilon0 = 1./(4. * np.pi)   #permittivity of free space
epsilon  = -1.               #sign of electron charge

#parameters
omega = 1.5 * 527646296.344         #wave frequency
B     = 30.                         #background B strength
n_e   = 1.                          #electron number density cm^-3
w_T   = 10.                         #dimensionless electron temp. k_B T / m c^2
theta = np.pi / 3.                  #observer angle

#derived quantities

k_perp = omega / c * np.sin(theta)                  # wavevector perp component n = 1 approximation
k_z    = omega / c * np.cos(theta)                  # wavevector parallel comp. n = 1 approximation
omega_p = np.sqrt(n_e * e**2. / (m * epsilon0))     # plasma frequency
omega_c = e * B / (m * c)                           # cyclotron frequency
lambd   = k_perp**2. * w_T**2. / (2. * omega_c**2.) # just a defined parameter


