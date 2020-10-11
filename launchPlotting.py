import numpy as np
import matplotlib.pyplot as plt
import electronVehicleOptimisation_engineComparison
import pickle

plt.close()

fig_vehicleTotalMass = plt.figure()
ax_vehicleTotalMass = fig_vehicleTotalMass.add_subplot(1, 1, 1)

fig_vehicleLength = plt.figure()
ax_vehicleLength = fig_vehicleLength.add_subplot(1, 1, 1)

fig_vehicleCost = plt.figure()
ax_vehicleCost = fig_vehicleCost.add_subplot(1, 1, 1)

fig_vehicleMassDry = plt.figure()
ax_vehicleMassDry = fig_vehicleMassDry.add_subplot(1, 1, 1)

thrustLabels = [
    'De-Laval',
    'Constant-design-thrust',
    'Perfect-nozzle',
    'Aerospike'
    ]

plotLabels = [
    'De Laval, nozzle design alt: stage 1 = 15e3m, stage 2 = 25e3m',
    'Constant thrust, nozzle design alt: stage 1 = 15e3m, stage 2 = 25e3m',
    'Perfect nozzle, nozzle always ideally expanded',
    'Aerospike, nozzle design alt: stage 1 = 25e3m, stage 2 = 25e3m'
]

for thrustLabel in thrustLabels:
    # Load and plot the optimisation
    pickleSaveName = "00_optiSols_{0}.pickle".format(thrustLabel)
    opti_sols = pickle.load(open(pickleSaveName, "rb"))
    print(opti_sols)


    def collectVariableFromOpti_sols(variable):
        variable = np.array([opti_sol[variable] for opti_sol in opti_sols])
        deltaV_setpoints = np.array([opti_sol['deltaV_setpoint'] for opti_sol in opti_sols])
        return deltaV_setpoints, variable


    deltaVs, vehicle_mass_totals = collectVariableFromOpti_sols('vehicle_mass_total')
    deltaVs, vehicle_mass_propellant = collectVariableFromOpti_sols('vehicle_mass_propellant')
    deltaVs, vehicle_mass_drys = collectVariableFromOpti_sols('vehicle_mass_dry')

    ax_vehicleMassDry.plot(deltaVs, vehicle_mass_drys, label=thrustLabel, marker='o')

    ax_vehicleTotalMass.plot(deltaVs, vehicle_mass_totals, label=thrustLabel, marker='o')

    vol = vehicle_mass_propellant/1000.0
    vol_vehicle = vol*1.2 # Say that the propellant tanks make up 80% of the volume.
    vehicleArea = np.pi*(1.2**2)/4
    length_vehicle = vol/vehicleArea
    ax_vehicleLength.plot(deltaVs, length_vehicle, label=thrustLabel, marker='o')

    fractionOfLength = length_vehicle/12.1
    normalCost = 5.0e6 # for a 12.1m vehicle
    vehicleCost = normalCost*fractionOfLength

    ax_vehicleCost.plot(deltaVs, length_vehicle, label=thrustLabel, marker='o')


ax_vehicleTotalMass.set_title('Electron launch vehicle, different engine types\n150kg payload, 8% structural ratio')
ax_vehicleTotalMass.set_ylabel('Total vehicle mass (kg)')
ax_vehicleTotalMass.set_xlabel('DeltaV setpoint (m/s')
ax_vehicleTotalMass.legend()
fig_vehicleTotalMass.savefig(fname='01_ElectronEngineComparison_vehicleMassTotal.png', dpi=200)

ax_vehicleLength.set_title('Electron launch vehicle, different engine types\n150kg payload, 8% structural ratio')
ax_vehicleLength.set_ylabel('Vehicle length (m)')
ax_vehicleLength.set_xlabel('DeltaV setpoint')
ax_vehicleLength.legend()
fig_vehicleLength.savefig(fname='01_ElectronEngineComparison_vehicleLength.png', dpi=200)



ax_vehicleCost.set_title('Electron launch vehicle, different engine types\n150kg payload, 8% structural ratio')
ax_vehicleCost.set_ylabel('Vehicle Cost (USD)')
ax_vehicleCost.set_xlabel('DeltaV setpoint')
ax_vehicleCost.legend()
fig_vehicleCost.savefig(fname='01_ElectronEngineComparison_vehicleCost.png', dpi=200)


ax_vehicleMassDry.set_title('Electron launch vehicle, different engine types\n150kg payload, 8% structural ratio')
ax_vehicleMassDry.set_ylabel('Vehicle mass dry (kg)')
ax_vehicleMassDry.set_xlabel('DeltaV setpoint')
ax_vehicleMassDry.legend()
fig_vehicleMassDry.savefig(fname='01_ElectronEngineComparison_vehicleMassDry.png', dpi=200)

