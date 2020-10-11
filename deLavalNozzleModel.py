# Written by Graham Bell to understand the basic performance of rocket nozzles and
# effects of various parameters no rocket flight performance
# 20/Jan/2019.

import numpy as np
import matplotlib.pyplot as plt
from atmosphericModels import atmosPress_earth, atmosPress_mars
from rocketShip_class import rocketShip
from scipy.optimize import fsolve

class rocketEngine_deLaval:
    def __init__(self, engine_mdot, P_1=10e5, designAlt=None, P_2=None, v_1=0.0, atmos_model=atmosPress_earth,
                 overExpSeparation=None, propellantProps=None, **kwargs):
        self.P_1 = P_1 # Pa, Combustion chamber press

        # self.designAlt = designAlt # m, nozzle design altitude, calculates level of expansion.
        self.atmos_model = atmos_model
        if P_2 is None:
            # Must input design altitude to back calculate the design exit pressure, P_2
            self.designAlt = designAlt
            _, self.P_2, _ = self.atmos_model(alt=self.designAlt)
        else:
            self.P_2 = P_2

        if designAlt is not None:
            pass # self.designAlt =
        # if designAlt is None:
        #     def optimisationFunc(test_alt, P_2):
        #         _, P_2_try, _ = self.atmos_model(alt=test_alt)
        #         return P_2 - P_2_try
        #     self.designAlt = fsolve(func=optimisationFunc, x0=1000.0, args=self.P_2)[0]
        # else:


        self.v_1 = v_1 # m/s, inital velocity of the propellants.
        self.m_dot = engine_mdot

        # Over expansion should be around 0.25.
        self.overExpSeparation = overExpSeparation

        if propellantProps is not None:
            self.insertPropellantProps(**propellantProps)
        # self.designNozzle()

    def insertPropellantProps(self, c_p=10.0, c_v=8.0, R_sp=287.0, T_1=1000.0, **kwargs):
        self.c_p = c_p
        self.c_v = c_v
        self.k = self.c_p/self.c_v #
        self.R_sp = R_sp # Specific molecular weight
        self.T_1 = T_1 # Combustion chamber temperature.

    def calc_v_2(self, k=None, R_sp=None, T_1=None, P_2=None, P_1=None, v_1=None):
        if k is None:
            k=self.k
        if R_sp is None:
            R_sp = self.R_sp
        if T_1 is None:
            T_1 = self.T_1
        if P_2 is None:
            P_2 = self.P_2
        if P_1 is None:
            P_1 = self.P_1
        if v_1 is None:
            v_1 = self.v_1
        v_2 = np.sqrt((2*k)/(k-1.0) * R_sp * T_1 * (1.0-(P_2/P_1)**((k-1.0)/k)) + v_1)
        return v_2

    def calc_vol_sp_2(self, vol_sp_1=None, P_1=None, P_2=None, k=None):
        if vol_sp_1 is None:
            vol_sp_1 = self.vol_sp_1
        if P_1 is None:
            P_1 = self.P_1
        if P_2 is None:
            P_2 = self.P_2
        if k is None:
            k = self.k
        vol_sp_2 = vol_sp_1*(P_1/P_2)**(1.0/k)
        return vol_sp_2

    def designNozzle(self, suppressPrint=True):
        gammaPlus1 = self.k + 1
        gammaMinus1 = self.k - 1
        # Throat conditions @ M=1
        self.P_t = self.P_1 * (1.0 + (gammaMinus1/2))**(-self.k/gammaMinus1)
        self.T_t = self.T_1 * (1/(1+gammaMinus1/2))
        self.A_t = self.m_dot/self.P_t * np.sqrt((self.R_sp*self.T_t)/self.k)
        # Speed of sound in the throat.
        self.a_t = np.sqrt(self.k * self.R_sp * self.T_t)
        self.v_t = 1.0*self.a_t

        self.Mach_2 = np.sqrt((2/gammaMinus1)*((self.P_1/self.P_2)**(gammaMinus1/self.k)-1.0))
        self.A_2 = (self.A_t/self.Mach_2)*((1+(gammaMinus1/2)*self.Mach_2**2)/(gammaPlus1/2))**(gammaPlus1/(2*gammaMinus1))
        self.T_2 = self.T_t*(1.0+gammaMinus1/2.0 * self.Mach_2**2)**-1
        self.a_2 = np.sqrt(self.k * self.R_sp * self.T_2)
        self.v_2 = self.Mach_2*self.a_2

        self.F_design = self.m_dot * self.v_2
        if suppressPrint is False:
            print('Mass flow rate of propellants:', self.m_dot)
            print('Chamber Press', self.P_1 / 1e5, 'bar')
            print()
            print('Throat Press', self.P_t/1e5, 'bar')
            print('Throat temp', self.T_t)
            print('Throat area:', self.A_t, ': Diameter:', 1000.0*np.sqrt(4.0*self.A_t/np.pi), 'mm')
            print()
            print('Exit temp', self.T_2)
            print('Exit speed of sound', self.a_2)
            print('Exit Area:', self.A_2)
            print('Area ratio:', self.A_2 / self.A_t)
            print('Nozzle exit velocity:', self.v_2)
            print('Exit Mach:', self.Mach_2)
            print()
            print('Nozzle Pressure Ratio', self.P_1/self.P_2)

    def areaRatioMachNum(self, Mach=1.2, gamma=1.4):
        gm1 = gamma-1
        gp1 = gamma+1
        a2_On_aStar = np.sqrt((1.0/ Mach**2.0)*(((2.0 + gm1* Mach**2.0) / gp1)**(gp1 / gm1)))
        return a2_On_aStar

    def offdesignNozzle(self, P_1_offDesign, P_3):
        gammaPlus1 = self.k + 1
        gammaMinus1 = self.k - 1

        # These properties are different.
        P_t_offDesign = P_1_offDesign * (1.0 + (gammaMinus1 / 2)*(1.0**2)) ** (-self.k / gammaMinus1)
        # These properties are the same.
        Mach_2_offDesign = self.Mach_2
        # There is a different mach number because the overall pressure ratio IS different.
        # Mach_2_offDesign = np.sqrt((2/gammaMinus1)*((P_1_offDesign/P_3)**(gammaMinus1/self.k)-1.0))

        # With reference to exit mach number and P_1 total condition.
        P_2_offDesign = P_1_offDesign * (1.0 + (gammaMinus1 / 2)*(Mach_2_offDesign**2)) ** (-self.k / gammaMinus1)

        v_2_offDesign = self.v_2
        engine_mdot_offDesign = self.A_t*P_t_offDesign/np.sqrt((self.R_sp*self.T_t)/self.k)

        if False:
            # Npw need to calculate P2. Have area.
            def findMachNumForExitArea(Mach, A_2_design):
                A_2_try = self.A_t * self.areaRatioMachNum(Mach=Mach, gamma=self.k)
                return A_2_try - A_2_design

            # Get the exit area mach number. Ideally expanded mach number.
            Mach_2 = fsolve(func=findMachNumForExitArea, x0=10.0, args=self.A_2)

        A_2_offDesign = np.zeros(shape=P_1_offDesign.shape)
        A_2_offDesign[:] = self.A_2

        Mach_2_offDesign = np.zeros(shape=P_1_offDesign.shape)
        Mach_2_offDesign[:] = self.Mach_2

        T_2_offDesign = np.zeros(shape=P_1_offDesign.shape)
        T_2_offDesign[:] = self.T_2

        a_2_offDesign = np.zeros(shape=P_1_offDesign.shape)
        a_2_offDesign[:] = self.a_2

        v_2_offDesign = np.zeros(shape=P_1_offDesign.shape)
        v_2_offDesign[:] = self.v_2

        # Detect the separation
        if self.overExpSeparation is not None:
            print('P_2/P_3', P_2_offDesign / P_3)
            # If P_2 is less than 30% of P_3.
            sepDetected = np.less_equal(P_2_offDesign, self.overExpSeparation*P_3)
            print('sepDetected:', sepDetected)
            # If there are any separated indicies, recalculate those properties.
            if np.any(sepDetected):
                P_2_offDesign[sepDetected] = P_3 * self.overExpSeparation
                Mach_2_offDesign[sepDetected] = np.sqrt((2 / gammaMinus1) * ((P_1_offDesign[sepDetected]/P_2_offDesign[sepDetected]) ** (gammaMinus1 / self.k) - 1.0))
                # print(MachTest)
                A_2_offDesign[sepDetected] = (self.A_t / Mach_2_offDesign[sepDetected]) * (
                            (1 + (gammaMinus1 / 2) * Mach_2_offDesign[sepDetected] ** 2) / (gammaPlus1 / 2)) ** (
                                       gammaPlus1 / (2 * gammaMinus1))
                # A_2_offDesign[A_2_offDesign>self.A_2] = self.A_2
                T_2_offDesign[sepDetected] = self.T_t * (1.0 + gammaMinus1 / 2.0 * Mach_2_offDesign[sepDetected] ** 2) ** -1
                a_2_offDesign[sepDetected] = np.sqrt(self.k * self.R_sp * T_2_offDesign[sepDetected])
                v_2_offDesign[sepDetected] = Mach_2_offDesign[sepDetected] * a_2_offDesign[sepDetected]


                # plt.plot(P_1_offDesign/1e5, P_2_offDesign/1e5)
                # plt.plot(P_1_offDesign/1e5, P_1_offDesign/P_2_offDesign)
                # plt.plot(P_1_offDesign/1e5, Mach_2_offDesign)
                # plt.plot(P_1_offDesign/1e5, A_2_offDesign)
                # plt.show()
                # exit(0)
        # New pressure thrust.
        F_press_offDesign = (P_2_offDesign - P_3) * A_2_offDesign
        F_mom_offDesign = engine_mdot_offDesign * v_2_offDesign

        return F_press_offDesign, F_mom_offDesign, sepDetected, engine_mdot_offDesign


    def exitDiameter(self, designAlts):
        D_2 = np.sqrt(4.0/np.pi*self.A_2)
        D_2 = [D_2]*len(designAlts)
        return D_2

    def throatDiameter(self, designAlts=None):
        D_t = np.sqrt(4.0/np.pi*self.A_t)
        if designAlts is not None:
            D_t = [D_t]*len(designAlts)
        return D_t

    def throatSurfaceArea(self, designAlt):
        D_t = self.throatDiameter()
        SA_t = np.pi*np.array(D_t)
        return SA_t
    #
    # def exitSurfaceArea(self):
    #     D_2 = self.exitDiameter(designAlts=designAlts)
    #     SA_2 = np.pi * D_2
    #     return SA_2

    def momentumThrust(self):
        self.F_mom = self.m_dot * self.v_2
        return self.F_mom

    def pressureThrust(self, alts, ignoreSeparation=False):
        P_3 = self.atmos_model(alts)[1]
        # Reset the design altitude pressure.
        if True:
            if self.overExpSeparation is not None and ignoreSeparation is False:
                # As array to take the first value and repeat that out.
                # This is because P_2 may be a single float from prior calculations, or pressureThrust may have already
                # been called. Thus it would be an array and numpy will try and duplicate it.
                try:
                    self.P_2 = np.asarray(self.P_2)[0]
                    self.v_2 = np.asarray(self.v_2)[0]
                    self.vol_sp_2 = np.asarray(self.vol_sp_2)[0]
                    self.A_2 = np.asarray(self.A_2)[0]
                except:
                    pass
                P_2 = np.repeat(self.P_2, len(P_3))
                v_2 = np.repeat(self.v_2, len(P_3))
                vol_sp_2 = np.repeat(self.vol_sp_2, len(P_3))
                A_2 = np.repeat(self.A_2, len(P_3))

                # Recalculate P_2, A_2, and v_2 for the separated condition.
                print(P_2/P_3)
                sepDetected = np.less(P_2/P_3, self.overExpSeparation)
                print('sepDetected:', sepDetected)
                P_2[sepDetected] = P_3[sepDetected]*self.overExpSeparation
                v_2[sepDetected] = self.calc_v_2(P_2=P_2[sepDetected])
                vol_sp_2[sepDetected] = self.calc_vol_sp_2(P_2=P_2[sepDetected])
                A_2 = self.m_dot*vol_sp_2/v_2

                # These guys are now functions of altitude.
                self.P_2 = P_2
                self.v_2 = v_2
                self.A_2 = A_2

            else:
                sepDetected = np.zeros(P_3.size, dtype=bool)
        self.F_press = (self.P_2 - P_3)*self.A_2
        self.sepDetected = sepDetected
        return self.F_press

    def totalThrust(self, altitudes=None):
        pressureThrust = self.pressureThrust(alts=altitudes)
        # Since the pressure thrust was evaluated first, it takes care of adjusting v_2 for the momentum thrust.
        momentumThrust = self.momentumThrust()
        self.F_total = pressureThrust + momentumThrust
        return self.F_total


if __name__ == '__main__':

    # Propellant Generation
    keroLox = {
        'T_1':1000, # combustion temp
        'c_p': 10, # Specific heat at const press.
        'c_v': 8, # Specific heat at const Vol.
        'R' : 287.0, # Molar mass of fuels.
        }

    Ch4Lox = {}

    # Rocket Generation
    engine1 = {
        'propellantProps': keroLox,
        'P_1': 20e5, # Pa, 20 bar chamber pressure
        'v_1': 0.0, # Initial velocity of the propellant gas
        'designAlt': 10e3, # m, nozzle exit design altitude
        'engine_mdot': 5.0,
        'Isp': 100.0 # Specific impulse of the engine.
    }

    # Basic Nozzle thurst profile in the atmosphere
    if False:
        engineObj = rocketEngine_deLaval(**engine1)
        # Design a nozzle
        engineObj.designNozzle()
        alts = np.linspace(0, 100e3, 200)
        thrust = engineObj.totalThrust(altitudes=alts)
        plt.plot(alts, thrust, label='total')
        plt.plot(alts, engineObj.pressureThrust(alts=alts), ls='--', label='pressure')
        plt.plot(alts, np.repeat(engineObj.momentumThrust(), len(alts)), ls='--', label='momentum')
        plt.xlabel('Alt (m)')
        plt.ylabel('Thrust (N)')
        plt.legend()
        plt.show()
        exit(0)

    # Test atmosphere Model
    if False:
        alts = np.linspace(0, 100e3, 200)
        trace = np.array([atmosPress_earth(alt=alt) for alt in alts])
        # T_kelvin, P_pascals, rho
        T_kelvin = trace[:,0]
        P_pascals = trace[:, 1]
        rho = trace[:, 2]

        plt.plot(alts, T_kelvin)
        plt.xlabel('Altitude (m)')
        plt.ylabel('Temperature (C)')
        plt.show()
        plt.close()

        plt.plot(alts, P_pascals)
        plt.xlabel('Altitude (m)')
        plt.ylabel('Pressure (Pa)')
        plt.show()
        plt.close()

        plt.plot(alts, rho)
        plt.xlabel('Altitude (m)')
        plt.ylabel('Density')
        plt.show()
        plt.close()

        exit(0)

    # Basic Nozzle thurst profile in the atmosphere
    if True:
        for separationType in np.linspace(0,0.8, 8):
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
            plt.plot(alts, thrust, label='total'+str(separationType))
            plt.plot(alts, engineObj.pressureThrust(alts=alts), ls='--', label='pressure'+str(separationType))
            plt.plot(alts, np.repeat(engineObj.momentumThrust(), len(alts)), ls='--', label='momentum'+str(separationType))
        plt.xlabel('Alt (m)')
        plt.ylabel('Thrust (N)')
        plt.legend()
        plt.show()
        exit(0)

    # Test Rocketship model
    if True:
        engine = rocketEngine_deLaval(**engine1)
        engine.designNozzle()
        rocket = rocketShip(engineObj=engine)
        # print(rocket.engine.totalThrust(altitudes=[7000]))
        # print(rocket.orbitalVelocity(g0=9.81, R_planet=6371e3, orbit_alt=400e3))

        rocketSettings = {
            'mass_dry': 100.0,
            'mass_fuel': 50.0,
            'm_dot': 5.0,
            'F_thrust_func': rocket.engine.totalThrust,
            'xSecArea': 0.5,
            'theta': np.pi/4.0,
            'ode_t_span': [0, 50.0],
            'x0': 0.0e3,
            'u0': 0.0,
            'maxStep': 0.05,
        }

        ode_sol_drag_gravity = rocket.atmosphericTrajectory(
            Cd_func=rocket.dummy_Cd_func,
            g=9.81, **rocketSettings)
        ode_sol_dragOnly = rocket.atmosphericTrajectory(
            Cd_func=rocket.dummy_Cd_func, **rocketSettings)
        ode_sol_gravityOnly = rocket.atmosphericTrajectory(
            g=9.81, **rocketSettings)
        ode_sol_freeSpace = rocket.atmosphericTrajectory(**rocketSettings)

        labels = ['Drag+Gravity', 'Drag only', 'Gravity only', 'Free space']
        solutions = [ode_sol_drag_gravity, ode_sol_dragOnly, ode_sol_gravityOnly, ode_sol_freeSpace]

        for odeSol, label in zip(solutions, labels):
            t = odeSol['t']
            x, u = odeSol['y'][0], odeSol['y'][1]
            altitude = x*np.sin(rocketSettings['theta'])
            plt.plot(altitude, u, label=label)
            plt.xlabel('Alt (m)')
            plt.ylabel('Velocity (m/s)')

        plt.legend()
        plt.show()
