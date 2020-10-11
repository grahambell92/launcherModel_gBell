Written by Graham Bell
20/Jan/2019.

# launcherModel
Code base for launch model classes for doing multi-stage rocket launch modelling. 

# Classes that do object modelling:
    deLavalNozzleModel:
        Models the deLaval nozzle performance
    perfectlyExpandedNozzleModel:
        Models the nozzle performance for a nozzle who's exit diameter is continouously
        adjusted to meet the local atmospheric conditions. This nozzle is the theoretical perfect nozzle.
        This class inherets methods from the deLavalNozzle.
    aerospikeNozzleModel:
        Models the aerospike nozzle, perfectly expanded nozzle until the design altitude then switches to the
        deLaval underexpanded thrust.

    rocketShip_class:
        General code structure for rocketship.
        Need to specify stage dicts:
        propellantProperties -> engine object -> stage1
        propellantProperties -> engine object -> stage2
        list of stage dicts are then sent into rocketShip:
        [stage1, stage2] -> rocketShip.

        The rocketShip class then is used to calculate the stage mass of the rocket, then the trajectory can be calculated.
        The trajectory modelling is currently at a constant angle. This is where the next stage of improvement is needed.

# Example scripts
deLavalVsAltitude_plots shows a complete example for extracting and plotting the different nozzle performances.
rocketLab_params shows a complete example for how to build a rocket and solve the trajectory.
stagingOptimiser_demo shows a complete example for how to numerically optimise a rocket stages to maximise/meet a delta-V.