# Written by Graham Bell.
# A class that returns the aerospike thrust vs altitude.
# Inheretts many methods and properties from the delaval nozzle class.

import numpy as np
import matplotlib.pyplot as plt
from deLavalNozzleModel import rocketEngine_deLaval
from scipy import optimize

class rocketEngine_aerospike(rocketEngine_deLaval):
    # Design a deLaval nozzle for > designAlt
    # Design a perfect nozzle for < designAlt

    def perfectNozzleThrust(self, altitudes):
        # Grab the altitude pressure
        P_3 = np.array([self.atmos_model(alt=alt)[1] for alt in altitudes])

        # Redesign the nozzle.
        self.P_2 = P_3
        self.designNozzle()

        # Pressure thrust is zero in a perfect nozzle. as P2 always=P3.
        pressureThrust = np.array(super().pressureThrust(alts=altitudes))
        # Since the pressure thrust was evaluated first, it takes care of adjusting v_2 for the momentum thrust.
        momentumThrust = np.array(super().momentumThrust())
        self.F_total = pressureThrust + momentumThrust
        return self.F_total

    def totalThrust(self, altitudes=None):
        thrusts = np.zeros(shape=len(altitudes))
        for index, altitude in enumerate(altitudes):
            # Return the thrust from a perfect nozzle below the design altitude (over-expanded region)
            if altitude < self.designAlt:
                thrusts[index] = self.perfectNozzleThrust(altitudes=[altitude])

            # Return the thrust from a normal delaval nozzle (under-expanded region)
            elif altitude > self.designAlt:
                _, self.P_2, _ = self.atmos_model(alt=self.designAlt)
                self.designNozzle()
                thrusts[index] = super().totalThrust(altitudes=[altitude])

        self.F_total = np.array(thrusts)
        return self.F_total

    def exitDiameter(self, designAlts=None):
        # Altitudes is not used, it depends on the design altitude.
        # D_4 is a type of D_t, but for the aerospike throat area.
        D_4 = np.sqrt(4.0/np.pi*(self.A_2))
        if designAlts is not None:
            D_4 = [D_4]*len(designAlts)
        return D_4

    def throatDiameter(self,designAlts=None):
        D_t = np.sqrt(4.0/np.pi*self.A_t)
        if designAlts is not None:
            D_t = [D_t]*len(designAlts)
        return D_t

    def throatGap(self):
        D_t = self.throatDiameter()
        D_4 = self.exitDiameter()
        return (D_4 - D_t)/2

    def throatSurfaceArea(self, designAlt=None):
        D_t = np.array(self.throatDiameter())
        D_4 = np.array(self.exitDiameter())
        SA_t = np.pi * (D_4 + D_t)
        return SA_t


    def calcDesignAltitude(self, desiredDiameter):
        def solveRHS(altitude):
            # Update the engine design
            self.designAlt = altitude
            _, self.P_2, _ = self.atmos_model(alt=self.designAlt)
            self.designNozzle()
            D_4 = self.exitDiameter(designAlts=[altitude])[0]
            error = D_4 - desiredDiameter
            return error

        designAlt = optimize.fsolve(func=solveRHS, x0=3e3)[0]
        return designAlt

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
        'P_1': 100e5, # Pa, 20 bar chamber pressure
        'v_1': 0.0, # Initial velocity of the propellant gas
        'designAlt': 50e3, # m, nozzle exit design altitude
        'engine_mdot': 7.8*9, # kg/s # Rocketlab rutherford mdot and num first stage engines.
    }
    rocket = rocketEngine_aerospike(**engine1)
    alts = np.linspace(100,50e3, 100)
    thrust = rocket.totalThrust(altitudes=alts)

    if False:
        plt.plot(thrust)
        plt.xlabel('Altitude (m)')
        plt.ylabel('Thrust (N)')
        plt.show()

    if True:
        # Plot the engine fully expanded diameter
        designAltitude = rocket.calcDesignAltitude(desiredDiameter=1.2)
        print('designAltitude', designAltitude)