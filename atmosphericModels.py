import numpy as np

def atmosPress_mars(alt=10.0, R=191.8):
    # Curves and layers taken from nasa.
    # https://www.grc.nasa.gov/www/k-12/airplane/atmosmrm.html
    if alt < 7000:
        # Troposphere
        T = -31.0 - 0.000998 * alt
        P = 0.699 * np.exp(-0.00009*alt)

    elif alt > 7000:
        # Lower stratosphere
        T = -23.4 - 0.00222 * alt
        P = 0.699 * np.exp(-0.00009 * alt)

    T_kelvin = T + 273.0
    P_pascals = P * 1000.0
    rho = P_pascals/(R*T_kelvin)
    out = {
        'T': T_kelvin,
        'P': P_pascals,
        'rho': rho
    }
    return T_kelvin, P_pascals, rho #out

def atmosPress_earth(alt=10.0, R=286.7):
    # Curves and layers taken from nasa.
    # asarray handles whether you pass this function a float or a list or an array.
    alt = np.asarray(alt)
    T = np.zeros(shape=alt.shape)
    P = np.zeros(shape=alt.shape)

    # T = np.asarray()
    # T, P = 0.0, 101.325

    # altitudes that are < 11e3
    alt_inds = np.less(alt, 11e3)
    T[alt_inds] = 15.04 - 0.00649*alt[alt_inds]
    P[alt_inds] = 101.29 * ((T[alt_inds] + 273.1) / 288.08) ** 5.256

    # altitudes that are >= 11e3 amd < 25e3
    alt_inds = np.greater_equal(alt, 11e3) & np.less(alt, 25e3)
    T[alt_inds] = -56.46
    P[alt_inds] = 22.65*np.exp(1.73-0.000157*alt[alt_inds])

    # altitudes that are >= 25e3
    alt_inds = np.greater_equal(alt, 25e3)
    T[alt_inds] = -131.21 + 0.00299 * alt[alt_inds]
    P[alt_inds] = 2.488 * ((T[alt_inds] + 273.1) / 216.6) ** -11.388

    T_kelvin = T + 273.0
    # P is in kPa
    P_pascals = P*1000.0
    rho = P_pascals/(R*T_kelvin)
    return T_kelvin, P_pascals, rho

def pressure_gascolumn(altitude, R=8.31432, lapseRate=-0.0065, P_sealevel=101325, T_sealevel=293.0):
    # lapse rate is temperature change per m. Units are k/m.
    g = 9.81 # m/s**2
    molarMass_air = 0.0289644 # kg/mol
    altitude_sealevel = 0.0
    P_altitude =  P_sealevel * (1 + lapseRate/T_sealevel*(altitude-altitude_sealevel))**((-g*molarMass_air)/(R*lapseRate))
    return P_altitude

def pressure_gascolumn_eng(altitude, R=8.31432, lapseRate=-0.0065, P_sealevel=101325, T_sealevel=293.0):
    # lapse rate is temperature change per m. Units are k/m.
    P_altitude = 101325 * (1-2.25577e-5 * altitude)**5.25588
    return P_altitude

if __name__ == '__main__':
    import matplotlib.pyplot as plt

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

    if True:
        altitudes = np.linspace(-5e3, 20e3, 400)
        plt.plot(altitudes, pressure_gascolumn(altitudes)/1000.0)
        plt.plot(altitudes, pressure_gascolumn_eng(altitudes)/1000.0)
        plt.plot(altitudes, atmosPress_earth(alt=altitudes)[1]/1000.0)
        print(atmosPress_earth(alt=5000))
        plt.grid()
        plt.show()