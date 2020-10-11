# Written by Graham Bell to understand the basic performance of rocket nozzles and effects of various parameters

import numpy as np
import matplotlib.pyplot as plt
from atmosphericModels import atmosPress_earth, atmosPress_mars
from rocketShip_class import rocketShip
from deLavalNozzleModel import rocketEngine_deLaval

class rocketEngine_perfectlyExpanded(rocketEngine_deLaval):

    def totalThrust(self, altitudes=None):
        # Grab the altitude pressure
        P_3 = np.array([self.atmos_model(alt=alt)[1] for alt in altitudes])

        # Redesign the nozzle.
        self.P_2 = P_3
        self.designNozzle()
        # Pressure thrust is zero in a perfect nozzle. as P2 always=P3.
        self.pressureThrust = super().pressureThrust(alts=altitudes, ignoreSeparation=True)
        # Since the pressure thrust was evaluated first, it takes care of adjusting v_2 for the momentum thrust.
        self.momentumThrust = super().momentumThrust()
        self.F_total = self.pressureThrust + self.momentumThrust
        return self.F_total

    def exitDiameter(self, designAlts):

        # Grab the altitude pressure
        P_3 = np.array([self.atmos_model(alt=alt)[1] for alt in designAlts])

        # Redesign the nozzle.
        self.P_2 = P_3
        self.designNozzle()
        D_2 = np.sqrt(4.0/np.pi*self.A_2)
        return D_2


    def throatDiameter(self, designAlt=None):
        D_t = np.sqrt(4.0/np.pi*self.A_t)
        return D_t



if __name__ == '__main__':

    # Propellant Generation
    keroLox = {
        'T_1': 1000, # combustion temp
        'c_p': 10, # Specific heat at const press.
        'c_v': 8, # Specific heat at const Vol.
        'R' : 287.0, # Molar mass of fuels.
        }

    # Rocket Generation
    engine1 = {
        'propellantProps': keroLox,
        'P_1': 20e5, # Pa, 20 bar chamber pressure
        'v_1': 0.0, # Initial velocity of the propellant gas
        'designAlt': 10e3, # m, nozzle exit design altitude
        'engine_mdot': 5.0, # kg/s
        'Isp': 100.0 # Specific impulse of the engine.
    }

    rocket = rocketEngine_perfectlyExpanded(**engine1)
    alts = np.linspace(100,50e3, 100)
    thrust = rocket.totalThrust(altitudes=alts)
    plt.plot(thrust)
    plt.ylabel('Thrust (N)')
    plt.xlabel('Altitude (m)')

    plt.show()
