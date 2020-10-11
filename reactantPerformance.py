# Calculate the reactant performance terms in the pressure to exhaust velocity rocket equation.
# Written by Graham Bell, 14/Jan/2019.

import numpy as np

def evaluateReactantPerformance(T_1=2000.0, M=20.0, gamma=1.4, density=1.4, **kwargs):
    R_universal = 8.314 # J K**-1 mol**-1
    if type(M) is list:
        # Average the list
        M = np.mean(M)
    R_sp = R_universal/M *1000.0 # Molar mass is in g/mol but needs to be kg.
    d = 2.0*gamma/(gamma-1.0)
    reactantPerformance = np.sqrt(d * R_sp * T_1)
    print('R_sp {0:.2f}, Performance {1:.2f}, Density performance {2:.2f}'.format(R_sp, reactantPerformance,
                                                                               reactantPerformance/density))


if __name__ == '__main__':
    gasConstR = 8.61

    air = {'Name': 'Air',
           'T_1':   273.16 + 15.0, # K, reservoir temperature/combustion temperature.
           'M':     28.97, # g/mol molar mass
           'gamma': 1.4, # ratio of specific heats.
           'density': 1.44, #kg/m**3
           }

    H202 = {'Name': 'H202',
            'T_1': 273.16 + 150.0,  # K, reservoir temperature/combustion temperature.
            'M': 34.01,  # g/mol molar mass
            'gamma': 1.241,  # ratio of specific heats.
            'density': 1.45*1000.0,  # kg/m**3
           }

    LOXH202 = {
        'Name': 'LOX/H202',
        'T_1': 273.16 + 15.0,  # K, reservoir temperature/combustion temperature.
        'M': [32.0, 34.01],  # g/mol molar mass
        'gamma': 1.4,  # ratio of specific heats.
        'density': 1.45*1000.0,  # kg
        }

    H202Kero = {
        'Name': 'H202/Kero',
        'T_1': 2975.0,  # K, reservoir temperature/combustion temperature.
        'M': [34.01, 167.5],  # g/mol molar mass
        'gamma': 1.31,  # ratio of specific heats.
        'density': 1.66*1000.0,  # kg
        }

    LOXKero = {
        'Name': 'LOX/Kero',
        'T_1': 3670.0,  # K, reservoir temperature/combustion temperature.
        'M': 167.5,  # g/mol molar mass
        'gamma': 1.22,  # ratio of specific heats.
        'density': np.mean([1.141, 0.86])*1000.0,  # kg
        }

    LOXEthanol = {
        'Name': 'LOX/Ethanol',
        'T_1': 3420.0,  # K, reservoir temperature/combustion temperature.
        'M': [32.0, 46.07],  # g/mol molar mass
        'gamma': 1.21,  # ratio of specific heats.
        'density': np.mean([1.141 ,0.98]) * 1000.0,  # kg
    }

    LOXCH4 = {
        'Name': 'LOXCH4',
        'T_1': 3533.0,  # K, reservoir temperature/combustion temperature.
        'M': [32.0, 16.04],  # g/mol molar mass
        'gamma': 1.26,  # ratio of specific heats.
        'density': np.mean([1.141, 0.656]) * 1000.0,  # kg
        }

    LOXH2 = {
        'Name': 'LOXH2',
        'T_1': 2959.0,  # K, reservoir temperature/combustion temperature.
        'M': [32.0, 4.0],  # g/mol molar mass
        'gamma': 1.26,  # ratio of specific heats.
        'density': np.mean([0.26]) * 1000.0,  # kg
        }

    reactants = [air, H202, LOXH202, LOXKero, LOXEthanol, LOXCH4, LOXH2]
    for reactant in reactants:
        print(reactant['Name'])
        evaluateReactantPerformance(**reactant)