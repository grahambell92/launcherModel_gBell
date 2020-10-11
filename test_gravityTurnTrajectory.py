import matplotlib.pyplot as plt
import numpy as np
# from deLavalPerformance.deLavalNozzleModel import rocketEngine
from rocketShip_class import rocketShip
from rocketLab_params import electronRocket_deLaval
from scipy import optimize
import pickle

rocketShip1 = rocketShip(**electronRocket_deLaval)

print()
print('#'*60)

# Do the optimisaton
rocketShip1.maxStep = 2.0

rocketShip1.updateStageMasses([8.0, 5.0])
ode_sols = rocketShip1.stagingManager()

for ode_sol in ode_sols:
    print(ode_sol.keys())
    theta = ode_sol['y'][2]
    x, u = ode_sol['y'][0], ode_sol['y'][1]
    print('velocity', u)
    print('theta', np.rad2deg(theta))
    exit(0)
    horiPos = ode_sol['horiPos']
    vertiPos = ode_sol['vertiPos']
    # print(horiPos/1e3)
    # print(vertiPos/1e3)
    # exit(0)
    plt.plot(horiPos/1000, vertiPos/1000, marker='o')


    # theta = ode_sol['y'][2]
    # time = ode_sol['t']
    # plt.plot(time, theta, marker='o')
plt.show()


print(ode_sols)
exit(0)


if False:
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



