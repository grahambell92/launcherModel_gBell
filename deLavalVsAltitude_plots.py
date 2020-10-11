import itertools
import numpy as np
from atmosphericModels import atmosPress_earth
from deLavalNozzleModel import rocketEngine_deLaval
from perfectlyExpandedNozzleModel import rocketEngine_perfectlyExpanded
from aerospikeNozzleModel import rocketEngine_aerospike

import matplotlib.pyplot as plt
import os


# Colors for plotting.
tableau10 = [(31, 119, 180),(255,127,14), (44, 160, 44), (214, 39, 40), (148, 103, 189),
             (140, 86, 75), (227, 119, 194), (127, 127, 127), (188, 189, 34), (23, 190, 207)]
for index, color in enumerate(tableau10):
    r, g, b = color
    tableau10[index] = (r / 255., g / 255., b / 255.)


# Propellant Generation
keroLox = {
    'T_1': 1000,  # combustion temp
    'c_p': 10,  # Specific heat at const press.
    'c_v': 8,  # Specific heat at const Vol.
    'R': 287.0,  # Molar mass of fuels.
}

# Rocket Generation
engine1 = {
    'propellantProps': keroLox,
    'P_1': 30e5,  # Pa, 20 bar chamber pressure
    'v_1': 0.0,  # Initial velocity of the propellant gas
    'designAlt': 60e3,  # m, nozzle exit design altitude
    'm_dot': 5.0,
    'overExpSeparation': None
}

saveFolder = './figures/deLavalPerformance/'
os.makedirs(saveFolder, exist_ok=True)

# exit velocity vs ratio of specific heats.
if False:
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1)

    P_1s = np.linspace(10e5, 300e5, 1000)
    ks = np.linspace(1.05, 1.4, 8)
    engineObj = rocketEngine_deLaval(**engine1)
    for k in ks:
        v_2 = engineObj.calc_v_2(k=k, R_sp=287.0, T_1=3000, P_2=1e5,
                                    P_1=P_1s, v_1=0.0)
        ax1.plot(P_1s/1e5, v_2, label='k:' + str(k))
    ax1.set_xlabel('Combustion chamber pressure, $P_1$ (bar)')
    ax1.set_ylabel('Exit velocity, $v_2$ (m/s)')
    ax1.legend()
    plt.show()
    plt.savefig(saveFolder+'00_exitVelvsGamma.png', dpi=200)
    plt.close()

# Single plot showing the components of a bell nozzle.
if False:
    # First show the components of thrust for a given nozzle.
    fig = plt.figure()
    ax1_thrust = fig.add_subplot(1, 1, 1)
    ax2_press = ax1_thrust.twinx()

    lineColors = itertools.cycle(tableau10)

    alts = np.linspace(start=0.0, stop=100e3, num=100)

    P_3 = np.array([atmosPress_earth(alt=alt)[1] for alt in alts])

    ax2_press.plot(alts, P_3/1000, label='Atmos. press.', color=next(lineColors))
    ax2_press.set_ylabel('Pressure (kPa)')
    ax2_press.legend(loc=0)
    designAlt = 10.0e3

    engine1['designAlt'] = designAlt

    engineObj = rocketEngine_deLaval(**engine1)
    # Design a nozzle
    engineObj.designNozzle()

    thrust = engineObj.totalThrust(altitudes=alts)
    F_mom = np.repeat(engineObj.momentumThrust(), len(alts))
    F_press = engineObj.pressureThrust(alts=alts)

    ax1_thrust.plot(alts, thrust, ls='-', label='Total thrust', color=next(lineColors))
    ax1_thrust.plot(alts, F_mom,label='Thrust, momentum component', ls='--', color=next(lineColors))
    ax1_thrust.plot(alts, F_press, label='Thrust, pressure component', ls='--', color=next(lineColors))

    ax1_thrust.set_xlabel('Altitude (m)')

    ax1_thrust.legend(loc=0)
    ax1_thrust.set_ylabel('Thrust (N)')
    ax1_thrust.set_xlabel('Altitude (m)')
    # plt.show()
    # exit(0)
    plt.savefig(saveFolder + '01_thrustComponents.png', dpi=200)
    plt.close()

 # Now a sweep across different design altitudes.
if False:

    fig = plt.figure()
    ax1_thrust = fig.add_subplot(1, 1, 1)
    ax2_press = ax1_thrust.twinx()

    lineColors = itertools.cycle(tableau10)

    alts = np.linspace(start=0.0, stop=100e3, num=100)

    P_3 = np.array([atmosPress_earth(alt=alt)[1] for alt in alts])

    ax2_press.plot(alts/1000, P_3/1000, ls='--', label='Atmos. press.')
    ax2_press.set_ylabel('Pressure (kPa)')

    designAlts = np.linspace(100, 20e3, 10)
    for designAlt in designAlts:
        engine1['designAlt'] = designAlt
        engineObj = rocketEngine_deLaval(**engine1)
        # Design a nozzle
        engineObj.designNozzle()
        thrust = engineObj.totalThrust(altitudes=alts)
        ax1_thrust.plot(alts/1000, thrust/1000, ls='-', label='Thrust, design alt: {:.1f}'.format(designAlt),
                        color=next(lineColors))

        # F_mom = np.repeat(engineObj.momentumThrust(), len(alts))
        # F_press = engineObj.pressureThrust(alts=alts)
        # ax1_thrust.plot(alts, F_mom,label='T_mom', ls='--')
        # ax1_thrust.plot(alts, F_press, label='T_press', ls='--')

    ax1_thrust.set_xlabel('Altitude (km)')
    ax1_thrust.set_ylabel('Thrust (kN)')

    ax2_press.legend(loc=0)

    ax1_thrust.legend(loc=0)
    # plt.show()
    # exit(0)
    plt.savefig(saveFolder+'02_designAltitude.png', dpi=200)
    plt.close()

# Bell nozzle compared to continuously changing exit area nozzle.
if False:
    fig = plt.figure()
    ax1_thrust = fig.add_subplot(1, 1, 1)
    ax2_press = ax1_thrust.twinx()

    lineColors = itertools.cycle(tableau10)

    alts = np.linspace(start=0.0, stop=100e3, num=100)

    P_3 = np.array([atmosPress_earth(alt=alt)[1] for alt in alts])

    ax2_press.plot(alts/1000, P_3/1000, ls='--', label='Atmos. press.')
    ax2_press.set_ylabel('Pressure (kPa)')

    # Design a bell fixed area nozzle
    designAlt = 10e3
    engine1['designAlt'] = designAlt
    engineObj = rocketEngine_deLaval(**engine1)
    engineObj.designNozzle()
    thrust = engineObj.totalThrust(altitudes=alts)
    ax1_thrust.plot(alts/1000, thrust/1000, ls='-', label='Thrust, design alt: {:.1f}'.format(designAlt),
                    color=next(lineColors))

    #### Now Design a continously varying exit area nozzle.
    engineObj = rocketEngine_perfectlyExpanded(**engine1)
    # Now vectorize P_2 to be equal to the atmospheric pressure over a range of altitudes
    thrust = engineObj.totalThrust(altitudes=alts)
    pressThrust = engineObj.pressureThrust
    momThrust = engineObj.momentumThrust
    ax1_thrust.plot(alts / 1000, thrust / 1000, ls='-', label='Continously A_2 varying nozzle'.format(designAlt), color=next(lineColors))
    ax1_thrust.plot(alts / 1000, pressThrust / 1000, ls='-', label='Continously A_2 varying nozzle_press'.format(designAlt),
                    color=next(lineColors))
    ax1_thrust.plot(alts / 1000, momThrust / 1000, ls='-', label='Continously A_2 varying nozzle_momentum'.format(designAlt),
                    color=next(lineColors))

    ax1_thrust.set_xlabel('Altitude (km)')
    ax1_thrust.set_ylabel('Thrust (kN)')

    ax2_press.legend(loc=0)

    ax1_thrust.legend(loc=0)
    # plt.show()
    # exit(0)
    plt.savefig(saveFolder+'03_designAltitudeVsContinouslyVaryNozzle.png', dpi=200)
    plt.close()

# Comparing Bell vs. Continously varying, vs. aerospike.
if True:
    # Bell nozzle compared to continuously changing exit area nozzle.
    dashes = ['-', '--', '-.', ':']
    dashCycler = itertools.cycle(dashes)

    fig = plt.figure()
    ax1_thrust = fig.add_subplot(1, 1, 1)

    alts = np.linspace(start=0.0, stop=100e3, num=100)

    # Plot the atmospheric pressure
    if True:
        ax2_press = ax1_thrust.twinx()
        lineColors = itertools.cycle(tableau10)
        P_3 = np.array([atmosPress_earth(alt=alt)[1] for alt in alts])
        ax2_press.plot(alts/1000, P_3/1000, ls='--', label='Atmos. press.')
        ax2_press.set_ylabel('Pressure (kPa)')

    # Design a bell fixed area nozzle
    # designAlt = 15e3
    # engine1['designAlt'] = designAlt

    engineModels = [rocketEngine_deLaval, rocketEngine_perfectlyExpanded, rocketEngine_aerospike]
    engineLabels = [ 'De Laval', 'Perfectly Expanded', 'Aerospike']

    if True:
        for engineModel, engineLabel in zip(engineModels, engineLabels):
            engineObj = engineModel(**engine1)
            thrust = engineObj.totalThrust(altitudes=alts)
            ax1_thrust.plot(alts/1000, thrust/1000, label=engineLabel,
                            color=next(lineColors),
                            ls=next(dashCycler))

        ax1_thrust.set_xlabel('Altitude (km)')
        ax1_thrust.set_ylabel('Thrust (kN)')

        ax2_press.legend(loc=0)

        ax1_thrust.legend(loc=0)
        plt.savefig(saveFolder+'05_thrustNozzleTypeComparison.png', dpi=200)
        plt.close()

    # Do an engine exit diameter plot based off the engineModels.
    if False:
        fig = plt.figure()
        ax1_thrust = fig.add_subplot(1, 1, 1)

        for engineModel, engineLabel in zip(engineModels, engineLabels):
            engineObj = engineModel(**engine1)
            #thrust = engineObj.totalThrust(altitudes=alts)
            engineDiam = engineObj.exitDiameter(altitudes=alts)
            ax1_thrust.plot(alts/1000, engineDiam, label=engineLabel,
                            color=next(lineColors),
                            ls=next(dashCycler))
        ax1_thrust.set_xlabel('Altitude (km)')
        ax1_thrust.set_ylabel('Engine Diameter (m)')

        ax2_press.legend(loc=0)

        ax1_thrust.legend(loc=0)
        # plt.show()
        fig.savefig(saveFolder + '06_NozzleDiameterTypeComparison.png', dpi=200)
        plt.close(fig)
        exit(0)

    # Throat surface/line area
    if True:
        fig = plt.figure()
        ax1 = fig.add_subplot(1, 1, 1)
        ax2 = ax1.twinx()
        ax3 = ax1.twiny()

        # Design altitude
        designAlts = np.linspace(start=0.0, stop=60e3, num=40)
        exitDiams = []
        for engineModel, engineLabel in zip(engineModels, engineLabels):
            lineColor = next(lineColors)
            for designAlt in designAlts:
                engine1['designAlt'] = designAlt
                engineObj = engineModel(**engine1)
                #thrust = engineObj.totalThrust(altitudes=alts)
                throat_SA = engineObj.throatSurfaceArea(designAlt=designAlt)
                ax1.semilogy(designAlt/1000, throat_SA, label=engineLabel,
                                color=lineColor,
                                marker='o')


                if engineLabel == engineLabels[2]:

                    throat_gap = engineObj.throatGap()
                    ax2.semilogy(designAlt / 1000, throat_gap, label=engineLabel,
                             color=lineColor,
                             marker='d')
                    exitDiams.append(engineObj.exitDiameter())
        ax1.set_ylim(1e-1, 1e1)
        ax1.set_xlabel('Design Altitude (km)')
        ax1.set_ylabel('Throat line area (m)')
        ax2.set_ylabel('Throat gap (m)')
        ax1_thrust.legend(loc=0)
        ax3.set_xlim(exitDiams[0], exitDiams[-1])
        ax3.set_label('Engine exit diameter (m)')

        plt.show()
        fig.savefig(saveFolder + '06_NozzleDiameterTypeComparison.png', dpi=200)
        plt.close(fig)
        exit(0)

exit(0)

# Single plot now showing the differences in non-separated flow and separated flow.
if True:
    # one for a few different design altitudes.
    designAlts = [3e3, 5e3, 10e3, 15e3, 20e3, 50e3]

    for designAlt in designAlts:
        fig = plt.figure()
        ax1_thrust = fig.add_subplot(1, 1, 1)
        ax2_press = ax1_thrust.twinx()

        lineColors = itertools.cycle(tableau10)


        alts = np.linspace(start=0.0, stop=100e3, num=100)
        P_3 = np.array([atmosPress_earth(alt=alt)[1] for alt in alts])

        ax2_press.plot(alts/1000.0, P_3/1000, color=next(lineColors), label='Atmos. press.')
        ax2_press.set_ylabel('Pressure (kPa)')

        sepNahs = [None]
        sepYehs = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
        sepYehorNahs = sepNahs + sepYehs
        sepYehBaseLabel = 'Seperation below $P_2/P_3 <$'
        sepYehLabels = [sepYehBaseLabel + str(sepYeh) for sepYeh in sepYehs]

        sepLabels = ['No separation model'] + sepYehLabels

        for sepYehorNah, sepLabel in zip(sepYehorNahs, sepLabels):
            engine1['designAlt'] = designAlt
            engine1['overExpSeparation'] = sepYehorNah
            engineObj = rocketEngine_deLaval(**engine1)
            # Design a nozzle
            engineObj.designNozzle()

            thrust = engineObj.totalThrust(altitudes=alts)
            # F_mom = np.repeat(engineObj.momentumThrust(), len(alts))
            # F_press = engineObj.pressureThrust(alts=alts)

            lineColor = next(lineColors)
            if sepYehorNah is not None:
                sepDetected = engineObj.sepDetected
                # Plot the separated as a dashed line.
                ax1_thrust.plot(alts[sepDetected]/1000.0, thrust[sepDetected]/1000.0, ls='--',
                                color=lineColor, label=sepLabel)
                # Plot the unseparated as a solid line.
                ax1_thrust.plot(alts[~sepDetected]/1000.0, thrust[~sepDetected]/1000.0, ls='-', color=lineColor)
                # ax1_thrust.plot(alts, F_mom,label='T_mom', ls='--')
                # ax1_thrust.plot(alts, F_press, label='T_press', ls='--')
            else:
                ax1_thrust.plot(alts/1000.0, thrust/1000.0, ls='-', label=sepLabel, color=lineColor)
        ax1_thrust.set_xlabel('Altitude (km)')
        ax1_thrust.set_xlim(0, 30e3/1000.0)
        ax1_thrust.set_ylabel('Thrust (kN)')
        ax1_thrust.set_ylim(bottom=0.0, top=12)

        ax2_press.legend(loc=8)
        ax1_thrust.legend(loc=5)
        plt.savefig(saveFolder+'04_separationPositionEffect_'+str(int(designAlt/1000))+'km' + '.png', dpi=200)
        plt.close()

        # plt.show()
