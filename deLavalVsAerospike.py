import numpy as np
from deLavalNozzleModel import rocketEngine_deLaval
from aerospikeNozzleModel import rocketEngine_aerospike
import matplotlib.pyplot as plt
import itertools

# Propellant Generation
keroLox = {
    'T_1': 1000,  # combustion temp
    'c_p': 10,  # Specific heat at const press.
    'c_v': 8,  # Specific heat at const Vol.
    'R': 287.0,  # Molar mass of fuels.
}
Molar_mass_ethanolLOX = 26e-3  # 29.7226e-3
R_star = 8.314459848
# From Joel
# m_dot_design = 1.33
# R = R_star / Molar_mass
# T_comb = 3500
# gamma = 1.2

ethanolLox = {
    'T_1': 3500,  # combustion temp
    'c_p': 10,  # Specific heat at const press.
    'c_v': 8,  # Specific heat at const Vol.
    'R': R_star/Molar_mass_ethanolLOX,  # Molar mass of fuels.
}

if False:



    engine1 = {
        'propellantProps': keroLox,
        'P_1': 66e5,  # Pa, 20 bar chamber pressure
        'v_1': 0.0,  # Initial velocity of the propellant gas
        # 'designAlt': 14e3, #-3231.0, #14.355e3,  # m, nozzle exit design altitude
        'P_2': 101325,
        'engine_mdot': 1.5,
        'Isp': 275.0,  # Specific impulse of the engine.
        'overExpSeparation': 0.25,
    }

    # engine1['overExpSeparation'] = separationRatio
    engineObj = rocketEngine_deLaval(**engine1)
    # Design a nozzle
    engineObj.designNozzle()
    alts = np.linspace(start=0, stop=100e3, num=200)
    thrust = engineObj.totalThrust(altitudes=alts)
    plt.plot(alts, thrust, label='total_delaval_sep')
    # plt.plot(alts, engineObj.pressureThrust(alts=alts), ls='--', label='pressure_' + label)
    # plt.plot(alts, np.repeat(engineObj.momentumThrust(), len(alts)), ls='--',
    #          label='momentum_' + label)
    plt.xlabel('Alt (m)')
    plt.ylabel('Thrust (N)')
    plt.legend()
    plt.show()
    exit(0)

if True:

    P_3 = 101325
    designNPR = 45.0
    P_2 = P_3
    P_1 = designNPR * P_3
    thrust_max = 3333.3
    Isp_max = 272.43
    engine_mdot_max = thrust_max/(Isp_max*9.81)

    engine1 = {
        'propellantProps': ethanolLox,
        'P_1': P_1,  # Pa, 20 bar chamber pressure
        'P_2': P_3,  # Engine exit pressure is equal to atmospheric pressure, P_3
        'v_1': 0.0,  # Initial velocity of the propellant gas
        # 'designAlt': -3231.0,  # -3231.0, #14.355e3,  # m, nozzle exit design altitude
        'engine_mdot': engine_mdot_max,
        'Isp': Isp_max,  # Specific impulse of the engine.
        'overExpSeparation': None,
    }
    engineObj = rocketEngine_deLaval(**engine1)
    # Design a nozzle to these specs.
    engineObj.designNozzle(suppressPrint=False)
    figsize = (15, 6)
    fig_thrust = plt.figure(figsize=figsize)
    ax1_press = fig_thrust.add_subplot(1,3,1)
    ax2_mom = fig_thrust.add_subplot(1,3,2)
    ax3_total = fig_thrust.add_subplot(1,3,3)

    fig_isp = plt.figure(figsize=figsize)
    ax4_press = fig_isp.add_subplot(1, 3, 1)
    ax5_mom = fig_isp.add_subplot(1, 3, 2)
    ax6_total = fig_isp.add_subplot(1, 3, 3)

    fig_thrustDiff = plt.figure(figsize=figsize)

    ax7_total_thrustDiff = fig_thrustDiff.add_subplot(1,1,1)

    if True:
        ax1_press.plot(engineObj.P_1 / 1e5, engineObj.pressureThrust(alts=[0]), marker='x', color='r', label='Design condition')
        ax2_mom.plot(engineObj.P_1/1e5, engineObj.momentumThrust(), marker='x', color='r', label='Design condition')
        ax3_total.plot(engineObj.P_1 / 1e5, engineObj.totalThrust(altitudes=[0]), marker='x', color='r', label='Design condition')

        ax4_press.plot(engineObj.P_1 / 1e5, engineObj.pressureThrust(alts=[0])/(9.81*engineObj.m_dot), marker='x', color='r', label='Design condition')
        ax5_mom.plot(engineObj.P_1 / 1e5, engineObj.momentumThrust() / (9.81 * engineObj.m_dot), marker='x',
                       color='r', label='Design condition')
        ax6_total.plot(engineObj.P_1 / 1e5, engineObj.totalThrust(altitudes=[0]) / (9.81 * engineObj.m_dot), marker='x', color='r', label='Design condition')

        ax7_total_thrustDiff.semilogy(engineObj.P_1 / 1e5, engineObj.totalThrust(altitudes=[0])-engineObj.totalThrust(altitudes=[0]), marker='x', color='r', label='Design condition')

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    lineColorCycler = itertools.cycle(colors)
    throttles = np.linspace(0.1, 1.5, num=1000)
    P_1s = P_1 * throttles

    engineObj.overExpSeparation = 0.0
    noSep_pressThrust, noSep_momThrust, sepDetected, noSep_engine_mdot_offdesign = engineObj.offdesignNozzle(P_1_offDesign=P_1s,
                                                                                           P_3=P_3, )
    noSep_totalThrust = noSep_pressThrust + noSep_momThrust
    for overSep in np.linspace(0.01, 1.0, num=10):
        engineObj.overExpSeparation = overSep
        # Now run the engine in all manner of other conditions.

        pressThrust, momThrust, sepDetected, engine_mdot_offdesign = engineObj.offdesignNozzle(P_1_offDesign=P_1s, P_3=P_3,) #


        totalThrust = pressThrust + momThrust
        lineColor = next(lineColorCycler)
        # print()
        # print('Momentum thrust')
        # print(momThrust)

        ## Pressure Thrust (loss
        # De Laval
        if False:
            ax1_press.plot(P_1s / 1e5, pressThrust,
                           label='De Laval', ls='-', color='b')

            # Aerospike
            ax1_press.plot(P_1s / 1e5, np.zeros(shape=P_1s.shape),
                           label='Aerospike', ls='-', color='r')

            ## Momentum Thrust
            # De Laval
            ax2_mom.plot(P_1s / 1e5, momThrust,
                         label='De Laval', ls='-', color='b')

            # Aerospike
            ax2_mom.plot(P_1s / 1e5, momThrust,
                         label='Aerospike', ls='-', color='r')

            ## Total Thrust
            # De Laval
            ax3_total.plot(P_1s / 1e5, totalThrust,
                           label='De Laval', ls='-', color='b')

            # Aerospike
            ax3_total.plot(P_1s / 1e5, momThrust,
                           label='Aerospike', ls='-', color='r')


        if True:
            ax1_press.plot(P_1s[~sepDetected] / 1e5, pressThrust[~sepDetected],
                           label='Separation factor: {0:.1f}'.format(overSep), ls='-', color = lineColor)
            ax1_press.plot(P_1s[sepDetected]/1e5, pressThrust[sepDetected], color=lineColor, ls='--',)

            ax2_mom.plot(P_1s[~sepDetected] / 1e5, momThrust[~sepDetected],
                         label='Separation factor: {0:.1f}'.format(overSep), ls='-', color=lineColor)
            ax2_mom.plot(P_1s[sepDetected] / 1e5, momThrust[sepDetected], color=lineColor, ls='--', )

            ax3_total.semilogy(P_1s[~sepDetected] / 1e5, totalThrust[~sepDetected],
                           label='Separation factor: {0:.1f}'.format(overSep), ls='-', color=lineColor)
            ax3_total.semilogy(P_1s[sepDetected] / 1e5, totalThrust[sepDetected], color=lineColor, ls='--', )

            ax4_press.plot(P_1s[~sepDetected] / 1e5, pressThrust[~sepDetected]/(9.81*engine_mdot_offdesign[~sepDetected]),
                           label='Separation factor: {0:.1f}'.format(overSep), ls='-', color = lineColor)
            ax4_press.plot(P_1s[sepDetected]/1e5, pressThrust[sepDetected]/(9.81*engine_mdot_offdesign[sepDetected]), color=lineColor, ls='--',)

            ax5_mom.plot(P_1s[~sepDetected] / 1e5, momThrust[~sepDetected]/(9.81*engine_mdot_offdesign[~sepDetected]),
                         label='Separation factor: {0:.1f}'.format(overSep), ls='-', color=lineColor)
            ax5_mom.plot(P_1s[sepDetected] / 1e5, momThrust[sepDetected]/(9.81*engine_mdot_offdesign[sepDetected]), color=lineColor, ls='--', )

            ax6_total.plot(P_1s[~sepDetected] / 1e5, totalThrust[~sepDetected]/(9.81*engine_mdot_offdesign[~sepDetected]),
                           label='Separation factor: {0:.1f}'.format(overSep), ls='-', color=lineColor)
            ax6_total.plot(P_1s[sepDetected] / 1e5, totalThrust[sepDetected]/(9.81*engine_mdot_offdesign[sepDetected]), color=lineColor, ls='--', )


            ax7_total_thrustDiff.semilogy(P_1s[~sepDetected] / 1e5, totalThrust[~sepDetected]-noSep_totalThrust[~sepDetected],
                           label='Separation factor: {0:.1f}'.format(overSep), ls='-', color=lineColor)
            ax7_total_thrustDiff.semilogy(P_1s[sepDetected] / 1e5, totalThrust[sepDetected]-noSep_totalThrust[sepDetected], color=lineColor, ls='--', )


        # plt.plot(P_1s / 1e5, momThrust)
        # plt.plot(P_1s / 1e5, momThrust+pressThrust, label='Separation factor: {0:.1f}'.format(overSep))
    ylabels = ['Pressure Thrust (N)', 'Momentum Thrust (N)', 'Total Thrust (N)']
    for ax, ylabel in zip([ax1_press, ax2_mom, ax3_total], ylabels):
        ax.grid()
        ax.set_xlabel('P_1 (bar)')
        ax.set_ylabel(ylabel)

    ax1_press.legend()

    ylabels = ['Pressure Isp (s)', 'Momentum Isp (s)', 'Total Isp (s)', 'Thrust difference (separated-non separated) (N)']
    for ax, ylabel in zip([ax4_press, ax5_mom, ax6_total, ax7_total_thrustDiff], ylabels):
        ax.grid()
        ax.set_xlabel('P_1 (bar)')
        ax.set_ylabel(ylabel)

    ax4_press.legend()
    ax7_total_thrustDiff.legend()

    fig_thrust.savefig('00_thrust_varyingChamberPressure_groundBased.png', dpi=200)
    fig_isp.savefig('01_isp_varyingChamberPressure_groundBased.png', dpi=200)
    fig_thrustDiff.savefig('02_deltaThrust_separatedToNonSep_varyingChamberPressure_groundBased.png', dpi=200)

    exit(0)

    plt.show()
    exit(0)

if True:
    engine1 = {
        'propellantProps': keroLox,
        'P_1': 66e5,  # Pa, 20 bar chamber pressure
        'v_1': 0.0,  # Initial velocity of the propellant gas
        'designAlt': 14e3, #-3231.0, #14.355e3,  # m, nozzle exit design altitude
        'engine_mdot': 1.5,
        'Isp': 275.0,  # Specific impulse of the engine.
        'overExpSeparation': None,
    }

    # NPR = chamber to atmosphere.
    # If you want to test altitudes, then you should vary the atmospheric pressure... but we can't
    # Instead you can vary the chamber pressure with a fixed atmospheric pressure
    # lower chamber pressure, lower altitude
    # higher chamber pressure, higher altitude.
    # So if we can achieve 66 bar chamber pressure, generating 10kN at sea level with 3 chambers.
    # Then a lower chamber pressure is equiv to a lower altitude.
    # at a 40 bar chamber pressure, the engine is at it's design condition. Over this it is under expanded. Under this it is over expanded.
    # At 10% throttle, that generates 6 bar chamber pressure. Which is strongly overexpanded.
    # 6 bar chamber pressure is an NPR=6/1=6
    # 66 bar chamber pressure is an NPR=66/1=66
    # The design alt is NPR=45. So engine chamber pressure at 45 bar.
    # So question is, what does the real engine look like that does NPR=6 at lift off, NPR = 45 at design alt.
    # Well liftoff at NPR=6 is just a 6 bar chamber pressure.
    # Design alt atmospheric pressure is then just NPR = P_chamb/P_atm = 45 = 6/P_atm = 6/45 = 0.1333 bar = 13,333 Pa = 14.355km
    # Design alt is then 45 = 66/p_atm = 66/45 = 1.46 bar. = -3231.66m

    # So the way to equate that with a real rocket is to do. Run it at x% of the max chamber pressure.
    # 6 chamber pressure @ 14km design alt.




    # rocket = rocketEngine_aerospike(**engine1)
    # alts = np.linspace(100, 50e3, 100)
    # thrust = rocket.totalThrust(altitudes=alts)



    for separationRatio in np.linspace(0, 0.8, 8):
        engine1['overExpSeparation'] = separationRatio
        engineObj = rocketEngine_deLaval(**engine1)
        # Design a nozzle
        engineObj.designNozzle()
        alts = np.linspace(start=0, stop=100e3, num=200)
        thrust = engineObj.totalThrust(altitudes=alts)
        plt.plot(alts, thrust, label='total_delaval_sep:' + str(separationRatio))
        # plt.plot(alts, engineObj.pressureThrust(alts=alts), ls='--', label='pressure_' + label)
        # plt.plot(alts, np.repeat(engineObj.momentumThrust(), len(alts)), ls='--',
        #          label='momentum_' + label)
    plt.xlabel('Alt (m)')
    plt.ylabel('Thrust (N)')
    plt.legend()
    plt.show()

    # if True:
    #     plt.plot(thrust)
    #     plt.xlabel('Altitude (m)')
    #     plt.ylabel('Thrust (N)')
    #     plt.show()
    #
    # if True:
    #     # Plot the engine fully expanded diameter
    #     designAltitude = rocket.calcDesignAltitude(desiredDiameter=1.2)
    #     print('designAltitude', designAltitude)
    exit(0)




# Basic Nozzle thurst profile in the atmosphere
if True:
    for separationType in np.linspace(0, 0.8, 8):
        engine2 = {
            'propellantProps': keroLox,
            'P_1': 20e5,  # Pa, 20 bar chamber pressure
            'v_1': 0.0,  # Initial velocity of the propellant gas
            'designAlt': 10e3,  # m, nozzle exit design altitude
            'engine_mdot': 5.0,
            'Isp': 100.0,  # Specific impulse of the engine.
            'overExpSeparation': separationType,
        }

        engineObj = rocketEngine_deLaval(**engine2)
        # Design a nozzle
        engineObj.designNozzle()
        alts = np.linspace(start=0, stop=100e3, num=200)
        thrust = engineObj.totalThrust(altitudes=alts)
        plt.plot(alts, thrust, label='total' + str(separationType))
        plt.plot(alts, engineObj.pressureThrust(alts=alts), ls='--', label='pressure' + str(separationType))
        plt.plot(alts, np.repeat(engineObj.momentumThrust(), len(alts)), ls='--',
                 label='momentum' + str(separationType))
    plt.xlabel('Alt (m)')
    plt.ylabel('Thrust (N)')
    plt.legend()
    plt.show()
    exit(0)