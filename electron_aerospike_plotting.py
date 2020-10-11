import numpy as np
import matplotlib.pyplot as plt
import pickle

fig_vehicleTotalMass = plt.figure()
ax_vehiclePayloadMass = fig_vehicleTotalMass.add_subplot(1, 1, 1)

fig_vehicleLength = plt.figure()
ax_vehicleLength = fig_vehicleLength.add_subplot(1, 1, 1)

fig_vehicleCost = plt.figure()
ax_vehicleCost = fig_vehicleCost.add_subplot(1, 1, 1)

fig_vehicleMassDry = plt.figure()
ax_vehicleMassDry = fig_vehicleMassDry.add_subplot(1, 1, 1)

# load electron-delaval and electron-aerospikeEnhanced
pickleOpenName = "00_optiSols_De-Laval.pickle"
opti_sols_deLaval = pickle.load(open(pickleOpenName, "rb"))
pickleOpenName = "01_optiSols_aerospike-payloadOptimisation.pickle"
opti_sols_aerospike_payloadOptimised = pickle.load(open(pickleOpenName, "rb"))


for optiSol_delaval, optiSol_aerospikeOptimised in zip(opti_sols_deLaval, opti_sols_aerospike_payloadOptimised):


    def collectVariableFromOpti_sols(opti_sols, variable):
        variable = np.array([opti_sol[variable] for opti_sol in opti_sols])
        deltaV_setpoints = np.array([opti_sol['deltaV_setpoint'] for opti_sol in opti_sols])
        return deltaV_setpoints, variable



    # deltaVs, vehicle_mass_totals = collectVariableFromOpti_sols('vehicle_mass_total')
    # deltaVs, vehicle_mass_propellant = collectVariableFromOpti_sols('vehicle_mass_propellant')
    deltaVs, vehicle_mass_drys_delaval = collectVariableFromOpti_sols(opti_sols_deLaval, 'vehicle_mass_dry')
    deltaVs, vehicle_mass_drys_aerospike = collectVariableFromOpti_sols(opti_sols_aerospike_payloadOptimised, 'vehicle_mass_dry')

    # deltaVs, vehicle_mass_payloads_delaval = collectVariableFromOpti_sols(opti_sols_deLaval, 'vehicle_mass_payload')
    deltaVs, vehicle_mass_payloads_aerospike = collectVariableFromOpti_sols(opti_sols_aerospike_payloadOptimised, 'vehicle_mass_payload')


    ax_vehicleMassDry.plot(deltaVs, vehicle_mass_drys_delaval, label='de Laval', marker='o')
    ax_vehicleMassDry.plot(deltaVs, vehicle_mass_drys_delaval, label='Aerospike', marker='o')

    # ax_vehiclePayloadMass.plot(deltaVs, vehicle_mass_payloads_delaval, label='de Laval', marker='o')
    ax_vehiclePayloadMass.plot(deltaVs, vehicle_mass_payloads_aerospike/150.0 * 100.0, label='Aerospike', marker='o')

    if False:
        vol = vehicle_mass_propellant/1000.0
        vol_vehicle = vol*1.2 # Say that the propellant tanks make up 80% of the volume.
        vehicleArea = np.pi*(1.2**2)/4
        length_vehicle = vol/vehicleArea
        ax_vehicleLength.plot(deltaVs, length_vehicle, label=thrustLabel, marker='o')

        fractionOfLength = length_vehicle/12.1
        normalCost = 5.0e6 # for a 12.1m vehicle
        vehicleCost = normalCost*fractionOfLength

        ax_vehicleCost.plot(deltaVs, length_vehicle, label=thrustLabel, marker='o')


# ax_vehiclePayloadMass.set_title('Electron launch vehicle, different engine types\n150kg payload, 8% structural ratio')
ax_vehiclePayloadMass.set_ylabel('Payload mass (kg)')
# ax_vehiclePayloadMass.hlines(150.0, xmin=8000, xmax=11e3)
ax_vehiclePayloadMass.set_xlabel('DeltaV setpoint (m/s')
ax_vehiclePayloadMass.legend()
fig_vehicleTotalMass.savefig(fname='02_ElectronEngineComparison_aerospikeVehiclePayloadOptimisation.png', dpi=200)

exit(0)

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

