import matplotlib.pyplot as plt
import numpy as np
# from deLavalPerformance.deLavalNozzleModel import rocketEngine
from rocketShip_class import rocketShip
from rocketLab_params import electronRocket_deLaval, electronRocket_perfectNozzle, \
    electronRocket_constThrust, electronRocket_aerospike
import itertools
from scipy import optimize
import pickle

dashes = ['-', '--', '-.', ':']
tableau10 = [(31, 119, 180), (255, 127, 14), (44, 160, 44), (214, 39, 40), (148, 103, 189),
                 (140, 86, 75), (227, 119, 194), (127, 127, 127), (188, 189, 34), (23, 190, 207)]
for index, color in enumerate(tableau10):
    r, g, b = color
    tableau10[index] = (r / 255., g / 255., b / 255.)

rocketParams = [
    # electronRocket_deLaval,
    # electronRocket_constThrust,
    electronRocket_perfectNozzle,
    # electronRocket_aerospike
]
thrustLabels = [
    # 'De-Laval',
    # 'Constant-design-thrust',
    'Perfect-nozzle',
    # 'Aerospike'
]

def main():

    # Optimisation of rocket stages
    print('#'*60)

    for rocketParam, thrustLabel in zip(rocketParams[:], thrustLabels[:]):
        rocketShip1 = rocketShip(**rocketParam)

        print()
        print('#'*60)
        print('Doing:', thrustLabel)

        pickleSaveName = "00_optiSols_{0}.pickle".format(thrustLabel)
        # Do the optimisaton
        rocketShip1.maxStep = 2.0
        minMR = 3.0
        maxMR = 10.0
        goal_deltaVs = np.linspace(8e3, 11.5e3, num=8)
        opti_sols = []
        initGuess = [4.0, 4.0]
        for goal_deltaV in goal_deltaVs:
            print('Solving deltaV:', goal_deltaV)
            deltaV_constraint = {"fun": rocketShip1.optimise_deltaVConstraint,
                                 "type": "eq",
                                 'args': (goal_deltaV,)
                                 }

            opti_sol = optimize.minimize(fun=rocketShip1.optimiseRocket_func,
                                         x0=initGuess,
                                         tol=1e-2,
                                         method='trust-constr',
                                         constraints=(deltaV_constraint),
                                         bounds = [(minMR, maxMR), (minMR, maxMR)],
                                         options={'gtol': 1.0,} # tolerence on meeing the deltaV constraint.
                                         # 1.0m/s within goal is sufficient.
                                         )
            # update the initial guess for the next optimisation at a different deltaV. Start closer to the goal.
            initGuess = opti_sol['x']

            rocketShip1.updateStageMasses(opti_sol['x'])
            opti_sol['deltaV_setpoint'] = goal_deltaV
            opti_sol['massRatios'] = opti_sol['x']

            opti_sol['vehicle_mass_total'] = rocketShip1.vehicle_mass_total
            opti_sol['mass_payload'] = rocketShip1.mass_payload
            opti_sol['vehicle_MR'] = rocketShip1.vehicle_MR
            opti_sol['vehicle_mass_dry'] = rocketShip1.vehicle_mass_dry
            opti_sol['vehicle_mass_propellant'] = rocketShip1.vehicle_mass_propellant
            opti_sol['vehicle_mass_total'] = rocketShip1.vehicle_mass_total
            opti_sol['vehicle_num_engines'] = rocketShip1.vehicle_num_engines
            for index, stage in enumerate(rocketShip1.stages):
                opti_sol['stage{0}_massRatio'.format(index)] = stage['massRatio']
                opti_sol['stage{0}_mass_dry'.format(index)] = stage['mass_dry']
                opti_sol['stage{0}_mass_fuel'.format(index)] = stage['mass_fuel']
                opti_sol['stage{0}_numEngines'.format(index)] = stage['numEngines']
            opti_sols.append(opti_sol)

            print('Mass ratio solution:', opti_sol['x'])
            print('Total vehicle mass:', rocketShip1.vehicle_mass_total)
            print()

        pickle.dump(opti_sols, open(pickleSaveName, "wb"))

        if False:
            print('Vehicle MR:', rocketShip1.vehicle_MR)
            print('Vehicle dry mass:', rocketShip1.vehicle_mass_dry)
            print('Vehcile propellant mass:', rocketShip1.vehicle_mass_propellant)
            print('Vehcile total mass:', rocketShip1.vehicle_mass_total)
            print('Vehicle engine number:', rocketShip1.vehicle_num_engines)


if __name__ == '__main__':
    main()

