from atmosphericModels import atmosPress_earth, atmosPress_mars
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

class rocketShip:
    def __init__(self, atmos_model_func=atmosPress_earth, stages=None, u0=0.0, x0=0.0, theta0=np.pi/2.0,
                 engineObj=None,
                 deltaV_total=None, Isp_generic=None, g0=9.81, mass_payload=10.0, maxStep=0.05,
                 **kwargs):

        self.atmos_model_func = atmos_model_func

        self.stages = stages
        if self.stages is not None:
            self.n_stages = len(stages)
        else:
            self.n_stages = None
        self.u0 = u0
        self.x0 = x0
        self.theta0 = theta0
        self.deltaV_total=deltaV_total
        self.Isp_generic = Isp_generic
        self.g0 = g0
        self.mass_payload = mass_payload

        self.maxStep = maxStep
        # self.maxStep = 10.0


        if engineObj is not None:
            self.engine = engineObj

    def setupEngine(self, engineObj):
        self.engine = engineObj


    def optimalStaging_unrestricted(self):

        # This method does the unrestricted staging so that the isp's don't need to be the same.
        # The method is based on the sections 11.7 in:
        # Orbital Mechanics for Eng Students, page 570.

        def massRatio(exhaust_vel, lagrangeMultiplier, structural_ratio):
            mr = (exhaust_vel*lagrangeMultiplier - 1.0)/ (exhaust_vel*structural_ratio*lagrangeMultiplier)
            return mr

        def optimiseRocketEq(lagrangeMultiplier):
            # Taken from page 576 of the book 'Orbital Mechanics for Eng students'.
            exhaust_vels=[]
            structRatios = []
            termA = 0.0
            termC = 0.0
            for stage in self.stages:
                exhaust_vel = stage['engineObj'].v_2/1000.0
                exhaust_vels.append(exhaust_vel)
                structRatio = stage['structural_ratio']
                structRatios.append(structRatio)

                termA += exhaust_vel*np.log(exhaust_vel*lagrangeMultiplier-1.0)
                termC += exhaust_vel*np.log(exhaust_vel*structRatio)

            sum_ev = np.sum(exhaust_vels)
            termB = sum_ev*np.log(lagrangeMultiplier)
            RHS = termA - termB - termC - self.deltaV_total/1000.0

            # Equation from book example: 11.5; Orbital Mechanics for Eng Students, page 577.
            # RHS = 3.924*np.log(3.924*lagrangeMultiplier - 1.0) + 3.434*np.log(3.434*lagrangeMultiplier - 1.0) + \
            #     2.943*np.log(2.943*lagrangeMultiplier - 1.0) - 10.30*np.log(lagrangeMultiplier) + 7.5089 - 10.0

            if np.any(np.isnan(RHS)):
                RHS = 10.0
            return RHS

        lagrangeMultiplier = fsolve(func=optimiseRocketEq, x0=1.0)[0]
        print('Lagrangian Muliplier for stages', lagrangeMultiplier)

        cumulativeMass = 0.0
        for stage in self.stages[::-1]:
            print('')
            exhaust_vel = stage['engineObj'].v_2
            structRatio = stage['structural_ratio']
            mr = massRatio(exhaust_vel=exhaust_vel/1000.0, lagrangeMultiplier=lagrangeMultiplier,
                           structural_ratio=structRatio)
            print('Stage mass ratio:', mr)
            stage['massRatio'] = mr
            sr = stage['structural_ratio']
            mass_stage = (mr - 1.0)/(1.0-mr*sr) * (self.mass_payload + cumulativeMass)
            mass_empty = sr*mass_stage
            mass_propellant = mass_stage - mass_empty


            print('mass_empty:', mass_empty, 'kg')
            print('mass_propellant:', mass_propellant, 'kg')

            stage['mass_empty'] = mass_empty
            stage['mass_fuel'] = mass_propellant
            stage['mass_dry'] = mass_stage

            # stage['mass_dry'] = stage['mass_empty'] + cumulativeMass
            # cumulativeMass += stage['mass_empty'] + stage['mass_fuel']
            total_mass = mass_stage + mass_propellant

            # # This is probably a dodgy way of upping the thrust to provide the acceleration.
            stageThrust0 = total_mass * stage['accel0']
            # # Because this will be the thrust when the stage has separated, not the design alt.
            stage['numEngines'] = stageThrust0 / stage['engineObj'].F_design
            print('Num engines:', stage['numEngines'])

    def updateStageMasses(self, massRatios):
        # Mass ratios are assumed to be [stage1, stage2, stage_n]
        cumulative_payload = self.mass_payload
        for stage, massRatio in zip(self.stages[::-1], massRatios[::-1]):
            structRatio = stage['structural_ratio']

            stage['massRatio'] = massRatio

            # Derived from the book Rocket Propulsion and Spaceflight Dynamics.

            payload_frac = (massRatio*structRatio - 1.0)/(massRatio*(structRatio-1.0))

            mass_observedPayload = cumulative_payload

            M_0 = (mass_observedPayload)/payload_frac

            propellant_frac = (massRatio - 1.0)/massRatio
            mass_propellant = M_0 * propellant_frac - mass_observedPayload

            mass_structure = M_0 - mass_observedPayload - mass_propellant

            mass_dry = mass_structure + mass_observedPayload

            cumulative_payload += mass_structure + mass_propellant

            # print('mass_dry:', mass_dry, 'kg')
            # print('mass_propellant:', mass_propellant, 'kg')

            stage['mass_dry'] = mass_dry
            stage['mass_fuel'] = mass_propellant
            total_mass = mass_dry + mass_propellant

            # # This is probably a dodgy way of upping the thrust to provide the acceleration.
            stageThrust0 = total_mass * stage['accel0']
            # # Because this will be the thrust when the stage has separated, not the design alt.
            stage['numEngines'] = stageThrust0 / stage['engineObj'].F_design

        cumulativeDryMass = 0.0
        cumulativePropellantMass = 0.0
        cumulativeEngines = 0.0
        for stage in self.stages:
            # print('Stage', stage['stage'], ':')
            cumulativeDryMass += stage['mass_dry']
            cumulativePropellantMass += stage['mass_fuel']
            cumulativeEngines += stage['numEngines']
            # print('Dry mass', stage['mass_dry'])
            # print('Fuel mass', stage['mass_fuel'])
            # print('Number of engines:', stage['numEngines'])
            # print('')

        vehicle_MR = (cumulativeDryMass + cumulativePropellantMass)/cumulativeDryMass
        self.vehicle_MR = vehicle_MR
        self.vehicle_mass_dry = cumulativeDryMass
        self.vehicle_mass_propellant = cumulativePropellantMass
        self.vehicle_mass_total = cumulativeDryMass + cumulativePropellantMass
        self.vehicle_num_engines = cumulativeEngines
        # print('Vehicle propellant mass', self.vehicle_mass_propellant)
        # print('Vehicle dry mass', self.vehicle_mass_dry)

    def updateStaging(self):
        # This method updates the stage mass properties based on the newly set mass ratios.
        # It goes through each of the stages in reverse and cumulates the mass to build the rocket.

        # Now add the distant payload masses to each stage in reverse order to get the correct payload mass on each stage.
        # E.g stage 1 sees a fuelled stage 2 as a payload, but stage 2 just sees the mass_payload as the payload.
        cumulativeMass = self.mass_payload
        for stage in self.stages[::-1]:
            stage['mass_dry'] = stage['mass_empty'] + cumulativeMass
            cumulativeMass += stage['mass_empty'] + stage['mass_fuel']

        print('')

        cumulativeEmptyMass = 0.0
        cumulativePropellantMass = 0.0
        cumulativeEngines = 0.0
        for stage in self.stages:
            print('Stage', stage['stage'], ':')
            cumulativeEmptyMass += stage['mass_empty']
            cumulativePropellantMass += stage['mass_fuel']
            cumulativeEngines += stage['numEngines']
            print('Dry mass', stage['mass_dry'])
            print('Fuel mass', stage['mass_fuel'])
            print('Number of engines:', stage['numEngines'])
            print('')

        # self.vehicle_MR = MR
        self.vehicle_mass_empty = cumulativeEmptyMass
        self.vehicle_mass_propellant = cumulativePropellantMass
        self.vehicle_num_engines = cumulativeEngines

        # print('Vehicle mass ratio:', MR)
        print('Total vehicle fuel:', cumulativePropellantMass)
        print('Total vehicle dry mass:', cumulativeEmptyMass)
        print('Total vehicle engines:', cumulativeEngines)

    def optimalStaging_restricted(self, structural_ratio=0.15):
        # This is the maths for doing the restricted staging problem.
        # The result is the mass ratios of each of the stages that results in an optimum delta v for given Isp
        # (or other effective performance parameter in the class inputs.

        # The restricted staging means that all stages have the same Isp, which is not that realistic.
        # The maths from this section come from
        # Rocket Propulsion and Spaceflight Dynamics,

        n_stages = len(self.stages)
        deltaV_total = self.deltaV_total

        Isp = self.Isp_generic
        MR = np.exp(deltaV_total / (n_stages * Isp * self.g0))
        generic_frac_payload = ((1 - MR * structural_ratio) / (MR * (1 - structural_ratio))) ** n_stages
        frac_payloads = [generic_frac_payload] * n_stages

        for index, (stage, payload_frac) in enumerate(zip(self.stages, frac_payloads)):
            stage['payload_frac'] = payload_frac
            stage['structural_ratio'] = structural_ratio

            SR = stage['structural_ratio']
            mass_payload = self.mass_payload
            index_a = 1.0 / self.n_stages
            index_b = ((self.n_stages - index) / self.n_stages)

            mass_empty = (((1.0 - payload_frac ** index_a) * SR) / (payload_frac ** index_b)) * mass_payload
            mass_propellant = (((1.0 - payload_frac ** index_a) * (1.0 - SR)) / (
                    payload_frac ** index_b)) * mass_payload

            print('mass_empty:', mass_empty, 'kg')
            print('mass_propellant:', mass_propellant, 'kg')

            stage['mass_empty'] = mass_empty
            stage['mass_fuel'] = mass_propellant
            total_mass = mass_empty + mass_propellant

            # # This is probably a dodgy way of upping the thrust to provide the acceleration.
            stageThrust0 = total_mass * stage['accel0']
            # # Because this will be the thrust when the stage has separated, not the design alt.
            stage['numEngines'] = stageThrust0 / stage['engineObj'].F_design

        self.updateStaging()

        self.MR = np.exp(self.deltaV_total / (self.n_stages * self.Isp_generic * self.g0))
        print('Vehicle mass ratio:', self.MR)

    def orbitalVelocity(self, g0, R_planet, orbit_alt):
        u_e = R_planet * np.sqrt(g0/(R_planet + orbit_alt))

        # Books formula is the same as this one.
        # u_e = np.sqrt(6.67408e-11 * 5.972e24/(R_planet + orbit_alt))
        return u_e

    def dummy_F_thrust_func(self, alt):
        return 10e3 # N

    def dummy_Cd_func(self, velocity):
        return 0.1

    def stagingManager(self):
        u0 = self.u0
        x0 = self.x0
        theta0 = self.theta0

        t0 = 0.0
        ode_sols = []
        for stage in self.stages:
            print('stage engines', stage['numEngines'])
            print('stage dry', stage['mass_dry'])
            # print('stage mdot', stage['engineObj'].m_dot)
            thisStage = {
                **stage,
                'u0': u0,
                'x0': x0,
                'theta0': theta0,
                'm_dot': stage['engineObj'].m_dot,
                'F_thrust_func': stage['engineObj'].totalThrust
            }

            # print('running on this stage:')
            # print(thisStage)
            # print()
            ode_sol = self.atmosphericTrajectory(**thisStage)

            t = ode_sol['t']
            x, u = ode_sol['y'][0], ode_sol['y'][1]
            # Gravity turn includes theta
            theta = ode_sol['y'][2]

            x0 = x[-1]
            u0 = u[-1]
            theta0 = theta[-1]
            # Add the remainding time to the previous solution.
            ode_sol['t'] += t0
            ode_sols.append(ode_sol)
            t0 = t[-1]
        return ode_sols



    def atmosphericTrajectory(self, mass_dry, mass_fuel, F_thrust_func, xSecArea=None, Cd_func=None, g0=None,
                              theta0=np.pi/2.0, ode_t_span=[0, 50000.0], u0=0.0, x0=0.0, m_dot=1.0,
                              numEngines=1.0, **kwargs):

        # This is the tracked event for fuel out.
        def fuelRemaining(t, y):
            residual_fuel = mass_fuel - (t * m_dot*numEngines)
            return residual_fuel

        # Signal to end the integration
        fuelRemaining.terminal = True

        appendedTimes = []
        altitudes = []
        drags = []
        gravities = []
        residualFuels = []
        thrusts = []
        horiPos = []
        vertPos = []
        horiSteps = []
        VertSteps = []

        self.x_last = 0.0
        # ODE Solving Func.
        def trajODE_RHS_straightline(t, y):
            x, u = y
            appendedTimes.append(t)

            altitude = x * np.sin(theta)
            altitudes.append(altitude)

            atmos_conds = self.atmos_model_func(alt=altitude)
            dragTerm = 0.0
            if Cd_func is not None:
                rho_alt = atmos_conds[2]
                dragTerm = Cd_func(u) * 0.5 * rho_alt * u**2 * xSecArea
                drags.append(dragTerm)

            gravityTerm = 0.0
            if g is not None:
                gravityTerm = g * np.sin(theta)
                gravities.append(gravityTerm)

            residual_fuel = fuelRemaining(t, y)
            residualFuels.append(residual_fuel)

            thrust = F_thrust_func([altitude])
            thrusts.append(thrust)

            # Straight acceleration
            dudt = (numEngines*thrust - gravityTerm - dragTerm)/(mass_dry + residual_fuel)
            dxdt = u

            return [dxdt, dudt]

        def trajODE_RHS_gravityTurn(t, y):
            x, u, beta = y
            appendedTimes.append(t)

            # theta=np.pi/2
            # altitude = x * np.sin(theta)
            # altitudes.append(altitude)

            x_step = x - self.x_last
            self.x_last = x

            horiz_step = x_step*np.cos(beta)
            vert_step = x_step*np.sin(beta)

            horiPos.append(np.sum(horiSteps) + horiz_step)
            vertPos.append(np.sum(horiSteps) + vert_step)

            horiSteps.append(horiz_step)
            VertSteps.append(vert_step)

            altitude = np.sum(VertSteps)
            altitudes.append(altitude)

            atmos_conds = self.atmos_model_func(alt=altitude)
            dragTerm = 0.0
            if Cd_func is not None:
                rho_alt = atmos_conds[2]
                dragTerm = Cd_func(u) * 0.5 * rho_alt * u ** 2 * xSecArea
                drags.append(dragTerm)

            gravityTerm = 0.0
            if g0 is not None:
                gravityTerm = g0 * np.sin(beta)
                gravities.append(gravityTerm)

            residual_fuel = fuelRemaining(t, y)
            residualFuels.append(residual_fuel)

            thrust = F_thrust_func([altitude])
            thrusts.append(thrust)

            # Gravity turn equation. Specify that the thrust is paralle to the velocity vector
            # theta is angle. Omega is angular rate
            gravityTerm = 9.81
            thrust = (numEngines * thrust - gravityTerm - dragTerm)
            currentMass = (mass_dry + residual_fuel)
            dudt =  (thrust / currentMass) - gravityTerm *np.cos(beta)
            dxdt = u
            dBetaDt = gravityTerm * np.sin(beta) / u
            return [dxdt, dudt, dBetaDt]

        # ode_sol = solve_ivp(fun=trajODE_RHS_straightline(), t_span=ode_t_span, y0=[x0, u0], max_step=self.maxStep, events=fuelRemaining)
        ode_sol = solve_ivp(fun=trajODE_RHS_gravityTurn, t_span=ode_t_span, y0=[x0, u0, theta0], max_step=self.maxStep, events=fuelRemaining)
        # Follow elliptical profile to orbit.
        


        # Put the extra properties into the ode_sol via interpolation
        if True:
            # Basic numpy interpolation to get back to ode_sol time steps.
            odeTime = ode_sol['t']
            appendedTimes_unique, uniqueInds = np.unique(appendedTimes, return_index=True)

            odeProps = [altitudes, drags, gravities, residualFuels, thrusts, horiPos, vertPos]
            odeProp_names = ['altitude', 'drag', 'gravity', 'residualFuel', 'engineThrust', 'horiPos', 'vertiPos']
            for odeProp, name in zip(odeProps, odeProp_names):
                # print(name)
                try:
                    # Todo Fix the engineThrust returning an array. Aka sort out the altitude inputs and outputs
                    odeProp = np.array(odeProp).squeeze()
                    odeProp_unique = odeProp[uniqueInds]
                    odeProp_interp = np.interp(x=odeTime, xp=appendedTimes_unique, fp=odeProp_unique)
                    ode_sol[name] = odeProp_interp
                except:
                    print('Failed to insert', name, 'into ode_sol')
        return ode_sol

    def plot_trajectory(self, ode_sol):
        t = ode_sol['t']
        x, u = ode_sol['y'][0], ode_sol['y'][1]
        plt.plot(t, u)
        #plt.plot(t, (mass_init- m_dot*t))
        # plt.plot(t, x)
        plt.show()

    def optimise_deltaVConstraint(self, params, goal_deltaV_total=3000.0):
        massRatios = params[:2]
        print(goal_deltaV_total)
        self.updateStageMasses(massRatios)
        ode_sols = self.stagingManager()
        final_sol = ode_sols[-1]
        x, u = final_sol['y'][0], final_sol['y'][1]
        final_u = u[-1]
        error = final_u - goal_deltaV_total
        print('DeltaV error:', error)
        return error

    def optimise_vehiclePropellantMass_constraint(self, params, goal_vehiclePropellantMass=4000):
        massRatios = params[:2]
        payloadMass_suggested = params[2]
        self.mass_payload = payloadMass_suggested
        self.updateStageMasses(massRatios)
        error = self.vehicle_mass_propellant - goal_vehiclePropellantMass
        print('Propellant error:', error)
        return error


    def optimiseRocket_func(self, massRatios):
        # This is the method called by the optimiser, it updates the new masses and runs the rocket solution and
        # returns error from desired deltaV

        # print('Trying mass ratios', massRatios)
        payloadMass_suggested = params[2]
        self.mass_payload = payloadMass_suggested
        self.updateStageMasses(massRatios)
        # print('MR:', massRatios, 'final u:', final_u)
        return self.vehicle_mass_total

    def optimisePayload_func(self, params):
        print('trying', params)
        massRatios = params[:2]
        payloadMass_suggested = params[2]

        # update payload mass
        self.mass_payload = payloadMass_suggested
        self.updateStageMasses(massRatios)

        return -self.mass_payload

if __name__ == '__main__':

    # Take a look at the demo at the bottom of rocketLab_params, that is how to use this class.
    pass