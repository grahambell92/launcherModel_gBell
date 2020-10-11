import numpy as np
from deLavalNozzleModel import rocketEngine_deLaval
from perfectlyExpandedNozzleModel import rocketEngine_perfectlyExpanded
from aerospikeNozzleModel import rocketEngine_aerospike
from atmosphericModels import atmosPress_earth

# Propellant Generation
keroLox = {
    'T_1': 4000,  # combustion temp
    'c_p': 10,  # Specific heat at const press.
    'c_v': 9.0,  # Specific heat at const Vol.
    'R': 287.0,  # Molar mass of fuels.
}

rutherFord_mdot = 7.8 # kg/s
# Rocket Generation
engineDict_stage1_general = {
    'propellantProps': keroLox,
    'P_1': 100e5,  # Pa, 20 bar chamber pressure
    'v_1': 0.0,  # Initial velocity of the propellant gas
    # 'designAlt': 15e3,  # m, nozzle exit design altitude
}


engineDict_stage2_general = {
    'propellantProps': keroLox,
    'P_1': 100e5,  # Pa, 20 bar chamber pressure
    'v_1': 0.0,  # Initial velocity of the propellant gas
}

engineDict_stage1_delaval = {
    **engineDict_stage1_general,
    'engine_mdot': rutherFord_mdot,
    'designAlt': 15e3,  # m, nozzle exit design altitude
    'overExpSeparation': 0.25
}

engineDict_stage2_delaval = {
    **engineDict_stage1_general,
    'engine_mdot': rutherFord_mdot, # kg/s
    'designAlt': 25e3,  # m, nozzle exit design altitude, # 25e3m produces an engine w/ exit diam of 1.2m
    'overExpSeparation': 0.25
}

engineDict_stage1_aerospike = {
    **engineDict_stage1_general,
    'engine_mdot': 9*rutherFord_mdot, # I want one engine with 9x the flow rate of equiv rocket lab vehicle.
    'designAlt': 25e3,  # m, nozzle exit design altitude
}

engineDict_stage2_aerospike = {
    **engineDict_stage2_general,
    'engine_mdot': 1*rutherFord_mdot, # kg/s
    'designAlt': 25e3,  # m, nozzle exit design altitude
}



###### Assign Engine Properties ########

engine1_obj_deLaval = rocketEngine_deLaval(**engineDict_stage1_delaval)
engine2_obj_deLaval = rocketEngine_deLaval(**engineDict_stage2_delaval)

engine1_obj_perfectNozzle = rocketEngine_perfectlyExpanded(**engineDict_stage1_delaval)
engine2_obj_perfectNozzle = rocketEngine_perfectlyExpanded(**engineDict_stage2_delaval)

engine1_obj_constThrust = rocketEngine_deLaval(**engineDict_stage1_delaval)
designThrust = engine1_obj_constThrust.totalThrust(altitudes=[engine1_obj_deLaval.designAlt])
def constThrust(altitudes):
    # Alt is not used.
    return designThrust

engine1_obj_constThrust.totalThrust = constThrust


engine2_obj_constThrust = rocketEngine_deLaval(**engineDict_stage2_delaval)
# Reassign the thrust function
designThrust = engine2_obj_constThrust.totalThrust(altitudes=[engine2_obj_constThrust.designAlt])
def constThrust(altitudes):
    # Alt is not used.
    return designThrust
engine2_obj_constThrust.totalThrust = constThrust


engine1_obj_aerospike = rocketEngine_aerospike(**engineDict_stage1_aerospike)
engine2_obj_aerospike = rocketEngine_aerospike(**engineDict_stage2_aerospike)

###### Stage Definitions ######

stage1 = {
    'stage': 1,
    # 'engineObj': None,
    'accel0': 1.4*9.81, # Initial acceleration of the stage.
    'theta0': np.pi/2.0,
    'ode_t_span': [0, 20e3],
    'structural_ratio': 0.08,
    'g0': 9.81,  # m/s**2

}

stage2 = {
    'stage': 2,
    # 'engineObj': None,
    'accel0': 2.0 * 9.81,  # Initial acceleration of the stage.
    'ode_t_span': [0, 20e3],
    # 'theta': np.pi/2.0,
    'structural_ratio': 0.08,
    'g0': 9.81,  # m/s**2

}

################## Stage properties assignment ##################

stage1_deLaval = {
    **stage1,
    'engineObj': engine1_obj_deLaval,
}

stage2_deLaval = {
    **stage2,
    'engineObj': engine2_obj_deLaval,
}

stage1_perfectThrust = {
    **stage1,
    'engineObj': engine1_obj_perfectNozzle,
}

stage2_perfectThrust = {
    **stage2,
    'engineObj': engine2_obj_perfectNozzle,
}

stage1_constThrust = {
    **stage1,
    'engineObj': engine1_obj_constThrust,
}

stage2_constThrust = {
    **stage2,
    'engineObj': engine2_obj_constThrust,
}

stage1_aerospike = {
    **stage1,
    'engineObj': engine1_obj_aerospike,
}

stage2_aerospike = {
    **stage2,
    'engineObj': engine2_obj_aerospike,
}


###### Stage Definitions ######

stages_deLaval = [stage1_deLaval, stage2_deLaval]
stages_perfectNozzle = [stage1_perfectThrust, stage2_perfectThrust]
stages_constThrust = [stage1_constThrust, stage2_constThrust]
stages_aerospike = [stage1_aerospike, stage2_aerospike]

electronRocket = {
    'u0': 50.0,
    'x0': 100.0, # 0.0,
    'theta0': np.deg2rad(0.001), # THeta is measured from the vertical, not horizontal.
    'atmos_model_func': atmosPress_earth,
    'mass_payload': 150.0,
    'Isp_generic': 311.0, # s # Need to update this to do it implicitly in the class. Sorry Dom.
}

electronRocket_deLaval = {
    **electronRocket,
    'stages': stages_deLaval
}

electronRocket_perfectNozzle = {
    **electronRocket,
    'stages': stages_perfectNozzle
}

electronRocket_constThrust = {
    **electronRocket,
    'stages': stages_constThrust
}

electronRocket_aerospike = {
    **electronRocket,
    'stages': stages_aerospike
}

if __name__ == '__main__':

    import matplotlib.pyplot as plt
    if True:
        rocketEngines = {
            'deLaval': engine1_obj_deLaval,
            'constant thrust': engine1_obj_constThrust,
            'Perfect nozzle': engine1_obj_perfectNozzle,
            'aerospike': engine1_obj_aerospike,

        }


        for engineName, engineObj in rocketEngines.items():
            print(engineName)
            alts = np.linspace(100, 50e3, 100)
            thrust = engineObj.totalThrust(altitudes=alts)
            plt.plot(alts, thrust, label=engineName)
        plt.xlabel('Altitude (m)')
        plt.ylabel('Thrust (N)')
        plt.legend()
        plt.show()


    if False:
        from rocketShip_class import rocketShip
        rocket1 = rocketShip(**electronRocket_deLaval)
        rocket1.deltaV_total = 9.0e3 # m/s
        rocket1.optimalStaging_restricted(structural_ratio=0.13)
        # rocket1.optimalStaging_unrestricted()

        exit(0)
        print('#'*60)

        rocket1 = rocketShip(**electronRocket_aerospike)
        rocket1.deltaV_total = 9000.0 # m/s
        rocket1.optimalStaging_restricted()
        rocket1.optimalStaging_unrestricted()