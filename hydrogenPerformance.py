import numpy as np
import matplotlib.pyplot as plt

LH2 = {
    'rho': 70.99, #kg/m**3
    'spec_energy': 119.3, # MJ/kg
    'energy_density': 8.491, # MJ/L
    'T_1': 123, # adiabatic flame temperature
    'c_p': 123,
    'c_v': 123,
    'R': 287.0,  # Molar mass of fuels.
}

LOX = {
    'rho': 1144.0, # kg/m**3
}

H2_LOX = {
    'rho': 123,
    'c_p': 123,
    'c_v': 123,
    'R': 287.0,  # Molar mass of fuels.
}

LNG = {
    'rho': 430.0, #kg/m**3
    'spec_energy': 53.6, # MJ/kg
    'energy_density': 22.2, # MJ/L
}

kero = {
    'rho': 800.0, #kg/m**3
    'spec_energy': 43.0, # MJ/kg
    'energy_density': 35.0, # MJ/L
}

# combustion modelling.
# mdot creates what pressure?
# mdot creates what temperature?
# Thrust v_2 = np.sqrt((2*k)/(k-1.0) * R_sp * T_1 * (1.0-(P_2/P_1)**((k-1.0)/k)) + v_1)

keroLox = {
    'T_1': 1000,  # combustion temp
    'c_p': 10,  # Specific heat at const press.
    'c_v': 8,  # Specific heat at const Vol.
    'R': 287.0,  # Molar mass of fuels.
}


# Rocket Generation
engine1 = {
    'propellantProps': keroLox,
    'P_1': 20e5,  # Pa, 20 bar chamber pressure
    'v_1': 0.0,  # Initial velocity of the propellant gas
    'designAlt': 10e3,  # m, nozzle exit design altitude
    'm_dot': 5.0,
    'Isp': 100.0  # Specific impulse of the engine.
}

# Basic Nozzle thurst profile in the atmosphere
if False:
    engineObj = rocketEngine_deLaval(**engine1)
    # Design a nozzle
    engineObj.designNozzle()
    alts = np.linspace(0, 100e3, 200)
    thrust = engineObj.totalThrust(altitudes=alts)

