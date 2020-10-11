import numpy as np

import matplotlib.pyplot as plt

if True:
    angles = np.linspace(start=0, stop=np.pi/2, num=50)
    radiuses_satellite = np.arange(100e3, 800e3, step=100e3)
    radius_earth = 6378.1e3 # radius of earth, m
    radiuses = radiuses_satellite + radius_earth
    G = 6.67E-11
    M1 = 5.972e24 # kg

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    for radius in radiuses:
        orbitalVelocity = np.sqrt(G*M1 / radius)
        deltaVs = 2*orbitalVelocity*np.sin(angles/2.0)
        ax.loglog(np.rad2deg(angles), deltaVs, label='orbit radius: {0} km'.format((radius - radius_earth)/1000.0))

    ax.set_ylabel('Required deltaV, m/s')
    ax.grid(True, which="both")
    ax.set_xlabel('Inclination change of circular orbit, deg')
    ax.legend()
    fig.savefig(fname='deltaVForInclinationChange.png', dpi=200)

if True:
    # Propellant needed for such deltaVs

    fig = plt.figure(figsize=(14, 6))
    ax_propMass = fig.add_subplot(1, 3, 1)
    ax_launchMass = fig.add_subplot(1, 3, 2)
    ax_launchCost = ax_launchMass.twinx()

    ax_costPerDegree = fig.add_subplot(1, 3, 3)

    spaceX_costPerkg = 5e3
    v_e = 20000.0 # m/s exit velocity for an ion thruster.
    m_payloads = np.arange(10.0, 280.0, step=30)
    m_propellants = np.linspace(1.0, 100.0, num=100)
    for m_payload in m_payloads:
        # m_struct = 0.08*m_propellants # 8% structural ratio per kg of propellant.
        m_struct = 80.0
        m_drymass = m_struct + m_payload
        m_init = m_drymass + m_propellants
        m_final = m_drymass
        MR = m_init/m_final
        deltaV = v_e * np.log(MR)
        # print(deltaV)
        ax_propMass.loglog(m_propellants, deltaV, label='Payload: {0:.0f} kg'.format(m_payload))
        ax_launchMass.loglog(deltaV, m_init, label='Payload: {0:.0f} kg'.format(m_payload))
        ax_launchCost.loglog(deltaV, m_init*spaceX_costPerkg, label='Payload: {0:.0f} kg'.format(m_payload))

        radius_earth = 6378.1e3  # radius of earth, m
        radius_satellite = 550e3 # Satellite altitude, m
        radius = radius_satellite + radius_earth
        G = 6.67E-11
        M1 = 5.972e24  # kg
        orbitalVelocity = np.sqrt(G * M1 / radius)
        deltaInclination = 2.0*np.arcsin(deltaV/(2.0*orbitalVelocity))
        costPerkg = (m_init*spaceX_costPerkg)/m_payload

        ax_costPerDegree.loglog(np.rad2deg(deltaInclination), costPerkg/1000.0, label='Payload: {0:.0f} kg'.format(m_payload))



    ax_propMass.grid(True, which="both")
    ax_propMass.set_ylabel('Delta V, m/s')
    ax_propMass.set_xlabel('Mass propellant, kg')
    ax_propMass.legend()

    ax_launchMass.grid(True, which="both", axis='x')
    ax_launchMass.set_xlabel('Delta V, m/s')
    ax_launchMass.set_ylabel('Mass propellant, kg')
    ax_launchMass.legend()

    ax_launchCost.grid(True, which="both")
    ax_launchCost.set_ylabel(r'Launch cost ($) at SpaceX rate ($5k/kg)')
    # ax_launchCost.set_xlabel('Mass propellant')

    ax_costPerDegree.grid(True, which="both")
    ax_costPerDegree.set_ylabel('Total Cost (k$)/kg inc. SpaceX base rate')
    ax_costPerDegree.set_yticks([10, 100])
    ax_costPerDegree.set_xlabel('Inclination Change, deg')
    ax_costPerDegree.legend()

    fig.suptitle('80kg space tug. Thruster exit velocity 20km/s, 550km orbit alt.')
    fig.tight_layout()
    fig.subplots_adjust(top=0.85)
    fig.savefig(fname='inclinationToLaunchCost.png', dpi=200)

