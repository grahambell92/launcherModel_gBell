
import numpy as np

costEstimations = {
    'Pressurizant Tank': (19.99465, 0.71253),
    'Fuel Tank': (19.99465, 0.71253),
    'Oxidiser Tank': (19.99465, 0.71253),
    'Thrust Cone': (2.79930, 0.91199),
    'Skirt': (2.79930, 0.91199),
    'Thermal Control': (2.79930, 0.91199),
    'Engine(s)': (31.48271, 0.78811),
    'Thrust Vector Control': (33.90978, 0.60977),
    # 'Pressurization System': (11.50618, 1.06948),
    'Pipes': (8.95877, 0.68815),
    'Valves': (8.95877, 0.68815),
    'Stage Harness': (27.45211, 0.44623),
    'Payload Adapter': (124.86209, 0.31031),
    'Payload Fairing': (4.09558, 0.96587),
    'Comms': (51.11253, 0.8),
    'Power': (56.13918, 0.66916),
    'Data Handling': (141.82428, 0.79249),
    'GNC': (69.05491, 0.82458),
    'Avionics Harness': (27.45211, 0.44623),
    # 'Attitude Control Module': (44.04074, 1.06207),
    'Interstage Structure': (6.70369, 0.68041),

}

if True:
    dryMass = 4000.0 # kg
    costDict = {}
    totalCost = 0.0
    for key, (a, b) in costEstimations.items():
        cost = a * (dryMass ** b) * 1e3 / 1e6  # into M of euros
        costDict[key] = cost
        totalCost += cost
    print('totalCost:', totalCost)
    for key, cost in costDict.items():
        print(key, cost/totalCost*100, '%')

exit(0)

if False:
    import matplotlib.pyplot as plt
    for dryMass in np.linspace(3000, 20000, num=10):
        costDict = {}
        totalCost = 0.0
        for key, (a, b) in costEstimations.items():
            cost = a*(dryMass**b)*1e3/1e6 # into M of euros
            costDict[key] = cost
            totalCost += cost
        plt.plot(dryMass, totalCost, marker='o')
    plt.show()
    exit(0)

