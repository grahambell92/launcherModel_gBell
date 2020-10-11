import numpy as np
from deLavalNozzleModel import rocketEngine_deLaval
from rocketShip_class import rocketShip
from atmosphericModels import atmosPress_earth

import matplotlib.pyplot as plt

# Propellant Generation
keroLox = {
    'T_1': 3800,  # combustion temp
    'c_p': 10,  # Specific heat at const press.
    'c_v': 9.5,  # Specific heat at const Vol.
    'R': 287.0,  # Molar mass of fuels.
    'k':1.2
}

# Rocket Generation
engine1 = {
    'P_1': 300e5,  # Pa, 20 bar chamber pressure
    'v_1': 0.0,  # Initial velocity of the propellant gas
    'designAlt': 1000,  # m, nozzle exit design altitude
    'm_dot': 5.0,
    'propellantProps': keroLox,

}

engine2 = {
    'P_1': 300e5,  # Pa, 20 bar chamber pressure
    'v_1': 0.0,  # Initial velocity of the propellant gas
    'designAlt': 1000,  # m, nozzle exit design altitude
    'm_dot': 5.0,
    'propellantProps': keroLox,

}


engine1 = rocketEngine_deLaval(**engine1)
engine1.designNozzle()
print('')

engine2 = rocketEngine_deLaval(**engine2)
engine2.designNozzle()
print('')

stage1 = {
    'stage': 1,
    'engineObj': engine1,
    'accel0': 1.0*9.81, # Initial acceleration of the stage.
    # 'F_thrust_func': constThrust,#engine1.totalThrust,
}

stage2 = {
    'stage': 2,
    # 'F_thrust_func': constThrust, #engine2.totalThrust,
}


singleStage = [stage1]
twoStage = [stage1, stage2]
threeStage = [stage1, stage2, stage2]
rocketConfigs = [singleStage, twoStage, threeStage]

for rocketStages in rocketConfigs:

    if False:
        # Working example 11.2
        n_stages = len(rocketStages)
        mass_payload = 10.0e3
        v_2 = 2000.0
        g0 = 9.81 # m/s**2
        Isp = 350.0

        # Pi_pl fraction
        frac_payload = 0.05

        # Structural ratio (eta)
        SR = 0.15
        # Mass ratio (n)
        # MR = (np.exp(deltaV_total/(n_stages*v_2)))**n_stages
        MR = 1/(frac_payload*(1.0-SR) + SR)

        # PR = (1.0-MR*SR)/(MR*(1-SR))

        delta_v = Isp*g0 * np.log(MR)
        m_0 = mass_payload/frac_payload

        mass_empty = SR*(m_0 - mass_payload)
        mass_propellant = m_0 - mass_empty - mass_payload

        print('Delta-V', delta_v, 'm/s')
        print('m_0:', m_0, 'kg')
        print('mass_empty:', mass_empty, 'kg')
        print('mass_propellant:', mass_propellant, 'kg')
    print('#')
    if True:
        # Working example 11.2
        n_stages = len(rocketStages)
        mass_payload = 10e3
        # v_2 = 4500.0
        # g0 = 9.81  # m/s**2
        # Isp = v_2/g0
        Isp = 350.0
        g0 = 9.81
        # delta_v = 7407.0 # m/s
        # Structural ratio (eta)
        SR = 0.15

        if True:
            # Solve by knowing frac_payload
            frac_payload = 0.05
            # Mass ratio (n)
            MR = 1/(frac_payload**(1/n_stages) * (1.0-SR) + SR)
            delta_v = Isp*g0 * np.log(MR**n_stages)
            print('Delta-V', delta_v)

        if False:
            # Solve by wanting desired delta_v
            delta_v = 8000.0 # m/s
            MR = np.exp(delta_v/(n_stages*Isp*g0))
            frac_payload = ((1-MR*SR)/(MR*(1-SR)))**n_stages

        for index, stage in enumerate(rocketStages):

            index_a = 1.0/n_stages
            index_b = ((n_stages-index)/n_stages)

            mass_empty = (((1.0-frac_payload**index_a)*SR)/(frac_payload**index_b)) * mass_payload
            mass_propellant = (((1.0-frac_payload**index_a)*(1.0-SR))/(frac_payload**index_b)) * mass_payload
            print('mass_empty:', mass_empty, 'kg')
            print('mass_propellant:', mass_propellant, 'kg')
