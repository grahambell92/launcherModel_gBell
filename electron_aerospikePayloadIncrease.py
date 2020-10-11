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

rocketShip1 = rocketShip(**electronRocket_aerospike)

# load de Laval information for vehicle masses and delta V setpoints.
pickleOpenName = "00_optiSols_De-Laval.pickle"
opti_sols_deLaval = pickle.load(open(pickleOpenName, "rb"))

pickleSaveName = "01_optiSols_aerospike-payloadOptimisation.pickle"


rocketShip1.maxStep = 2.0
initPayload = 150.0 # kg

initGuess = (4.0, 4.0, initPayload)
minMR = 3.0
maxMR = 10.0

opti_sols = []

for optiSol_deLaval in opti_sols_deLaval:

    deltaVSetpoint = optiSol_deLaval['deltaV_setpoint']
    vehiclePropellentMass = optiSol_deLaval['vehicle_mass_propellant']

    deltaV_constraint = {"fun": rocketShip1.optimise_deltaVConstraint,
                         "type": "eq",
                         'args': (deltaVSetpoint,) # arguments must have trailing comma.
                         }

    vehiclePropellantMass_constraint = {"fun": rocketShip1.optimise_vehiclePropellantMass_constraint,
                                        "type": "eq",
                                        'args': (vehiclePropellentMass,)
                         }

    # Given the same amount of propellant, how much additional payload can the vehicle carry?

    opti_sol = optimize.minimize(fun=rocketShip1.optimisePayload_func,
                                 x0=initGuess,
                                 tol=1e-2,
                                 method='trust-constr',
                                 constraints=(deltaV_constraint, vehiclePropellantMass_constraint),
                                 bounds = [(minMR, maxMR), (minMR, maxMR), (0, 300.0)],
                                 # options={'gtol': 1.0,} # tolerence on meeing the deltaV constraint.
                                 # 1.0m/s within goal is sufficient.
                                 )

    # update the initial guess for the next optimisation at a different deltaV. Start closer to the goal.
    initGuess = opti_sol['x']

    rocketShip1.updateStageMasses(opti_sol['x'])
    opti_sol['deltaV_setpoint'] = deltaVSetpoint
    opti_sol['massRatios'] = opti_sol['x']

    opti_sol['vehicle_mass_total'] = rocketShip1.vehicle_mass_total
    opti_sol['mass_payload'] = rocketShip1.mass_payload
    opti_sol['vehicle_MR'] = rocketShip1.vehicle_MR
    opti_sol['vehicle_mass_dry'] = rocketShip1.vehicle_mass_dry
    opti_sol['vehicle_mass_propellant'] = rocketShip1.vehicle_mass_propellant
    opti_sol['vehicle_mass_total'] = rocketShip1.vehicle_mass_total
    opti_sol['vehicle_num_engines'] = rocketShip1.vehicle_num_engines
    opti_sol['vehicle_mass_payload'] = rocketShip1.mass_payload

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


