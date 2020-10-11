import matplotlib.pyplot as plt
import numpy as np
# from deLavalPerformance.deLavalNozzleModel import rocketEngine
from rocketShip_class import rocketShip
# from deLavalPerformance.atmosphericModels import atmosPress_earth
from rocketLab_params import electronRocket_deLaval, electronRocket_perfectNozzle, \
    electronRocket_constThrust, electronRocket_aerospike
from perfectlyExpandedNozzleModel import rocketEngine_perfectlyExpanded
import itertools
from scipy import optimize

dashes = ['-', '--', '-.', ':']
tableau10 = [(31, 119, 180), (255, 127, 14), (44, 160, 44), (214, 39, 40), (148, 103, 189),
                 (140, 86, 75), (227, 119, 194), (127, 127, 127), (188, 189, 34), (23, 190, 207)]
for index, color in enumerate(tableau10):
    r, g, b = color
    tableau10[index] = (r / 255., g / 255., b / 255.)


if True:
    # Optimisation of rocket stages
    print('#'*60)

    colorCycler = itertools.cycle(tableau10)
    dashCycler = itertools.cycle(dashes)

    fig = plt.figure()
    plotCols = 3
    axThrust = fig.add_subplot(1, plotCols, 1)
    axTime = fig.add_subplot(1, plotCols, 2)
    axAlt = fig.add_subplot(1, plotCols, 3)



    rocketParams = [
        electronRocket_deLaval,
        electronRocket_constThrust,
        # electronRocket_perfectNozzle,
        electronRocket_aerospike
    ]
    thrustLabels = [
        'De Laval',
        'Constant design thrust',
        # 'Perfect nozzle',
        'Aerospike']
    for rocketParam, thrustLabel in zip(rocketParams[:], thrustLabels[:]):
        rocketShip1 = rocketShip(**rocketParam)
        lineStyle = next(dashCycler)

        print('#'*60)
        print('Doing:', thrustLabel)


        goal_deltaV = 9.5e3

        # rocketShip1.optimiseRocket_func(massRatios=[9.0,9.0])
        # exit(0)
        rocketShip1.maxStep = 10.0
        minMR = 2.0
        maxMR = 8.0
        opti_sol = optimize.minimize(fun=rocketShip1.optimiseRocket_func,
                                x0=[4.0, 4.0],
                                tol=1e-2,
                                args=(goal_deltaV),
                                method='Nelder-Mead'
                                # method = 'BFGS',
                                # method='SLSQP', bounds=((minMR, maxMR), (minMR, maxMR))
                                # method = 'TNC', bounds = ((minMR, maxMR), (minMR, maxMR))
                                # method = 'L-BFGS-B', bounds=((minMR, maxMR), (minMR, maxMR))
                                # method='CG',
                                )
        for row in opti_sol.items():
            print(row)

        if True:
            print('Vehicle MR:', rocketShip1.vehicle_MR)
            print('Vehicle dry mass:', rocketShip1.vehicle_mass_dry)
            print('Vehcile propellant mass:', rocketShip1.vehicle_mass_propellant)
            print('Vehcile total mass:', rocketShip1.vehicle_mass_total)
            print('Vehicle engine number:', rocketShip1.vehicle_num_engines)

        if True:
            massRatioSol = opti_sol['x']
            # Rerun the solution to get the ODE sols:
            rocketShip1.updateStageMasses(massRatios=massRatioSol)
            ode_sols = rocketShip1.stagingManager()

            for index, ode_sol in enumerate(ode_sols):
                t = ode_sol['t']
                x, u = ode_sol['y'][0], ode_sol['y'][1]
                alt = ode_sol['altitude']
                engineThrust = ode_sol['engineThrust']
                axThrust.plot(alt/1000.0, engineThrust, ls=lineStyle, label=thrustLabel, lw=2.0)
                axTime.plot(t, u, ls=lineStyle, label=thrustLabel, lw=2.0)
                axAlt.plot(t, alt / 1000.0, ls=lineStyle, label=thrustLabel, lw=2.0)

    axThrust.set_xlabel('alt (km)')
    axTime.set_xlabel('time (s)')
    axAlt.set_xlabel('time (s)')

    axThrust.set_ylabel('Thrust (N)')
    axTime.set_ylabel('Velocity (m/s)')
    axAlt.set_ylabel('Altitude (km)')

    axThrust.legend()
    axAlt.legend()
    plt.show()
    exit(0)