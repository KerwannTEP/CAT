##################################################
# Useful constants
##################################################

using SpecialFunctions

# Numerical Constants
# Calculation done with with Julia

const PI = 3.1415926535897

##################################################
# Constants of the problem
##################################################

const nbGlobularCluster = 10^5 # Simulation NBody

# Heggie et Hut (2003)
const parameterCoulomb = 0.11

# Coulomb logarithm
const logCoulomb = log(parameterCoulomb*nbGlobularCluster)

# Physical units
const _G = 1.0
const _M = 1.0
const _b = 1.0

# Dimensional units (Plummer potential)
const _v0 = sqrt(_G*_M/_b) # Velocity
const _E0 = -_G*_M/_b # Energy
const _L0 = sqrt(_G*_M*_b) # Action
const _F0 = (_G*_M*_b)^(-3/2) # Distribution function
const _Omega0 = sqrt(_G*_M/_b^3) # Frequency
const _rho0 = _M/(4*PI*_b^3/3) # Density


# Constant used to compute H(x,q) and the DF

const GAMMA_ca = gamma(9/2 - qCalc)
const GAMMA_db = gamma(1 - qCalc/2)
const GAMMA_bc = gamma(9/2 - qCalc/2)
const GAMMA_6q = gamma(6-qCalc)
