# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 18:56:05 2018

@author: JL
"""

import numpy as np
import matplotlib.pyplot as plt
import csv


def area_ratio_calc(P1, P2, gamma, T1):
    M = Mach_number(P1, P2, gamma)

    T = T1 * (1 + (gamma - 1) / 2 * M ** 2) ** (-1)
    A_R = ((gamma + 1) / 2) ** (-(gamma + 1) / (2 * (gamma - 1))) * (
                (1 + (gamma - 1) / 2 * M ** 2) ** ((gamma + 1) / (2 * (gamma - 1))) / M)

    return A_R, M, T


def exit_area(A_ratio, d_throat):
    A_throat = (d_throat / 2) ** 2 * np.pi

    A_exit = A_throat * A_ratio

    return A_exit


def exhaust_velocity(gamma, R, T, Mj):
    v_e = Mj * np.sqrt(gamma * R * T)

    return v_e


def Mach_number(P1, P2, gamma):
    M = np.sqrt(((P1 / P2) ** ((gamma - 1) / gamma) - 1) * 2 / (gamma - 1))

    return M


def thrust_output(P_variation, P_comb, gamma, P_design, thrust_design, A_exit, aerospike=False):
    if aerospike == True:
        thrust = thrust_design
    else:
        P_difference = P_design - P_variation
        thrust = thrust_design + P_difference * A_exit
        print(P_difference)

    return thrust


def choked_conditions(P1, T1, R, gam, m_dot):
    npr_crit = 1 / (2 / (gam + 1)) ** (gam / (gam - 1))

    rho_up = P1 / (R * T1)

    T_choked = T1 * (1 / npr_crit) ** ((gam - 1) / gam)
    c_choked = np.sqrt(R * gam * T_choked)
    rho_choked = rho_up * (1 / npr_crit) ** (1 / gam)

    A = m_dot / (rho_choked * c_choked)

    d_throat = np.sqrt(A / np.pi) * 2

    return T_choked, rho_choked, c_choked, d_throat


def set_constants(condition):
    R_star = 8.314459848
    P_SL = 101325
    T_SL = 288
    g_0 = 9.80665
    m_dot_design = 1.33

    alt_steps = 1000
    max_alt = 40e3
    altitudes = np.linspace(0, max_alt, alt_steps)
    #    altitudes = np.array([0,8.45e3,11.8e3,14.6e3,16.2e3])

    if condition == 'cold':
        Molar_mass = 28.9645e-3
        R = R_star / Molar_mass
        T_comb = 300
        P_comb = 500000
        gamma = 1.4
    elif condition == 'hot':
        Molar_mass = 26e-3  # 29.7226e-3
        R = R_star / Molar_mass
        T_comb = 3500
        P_comb = 6e6
        gamma = 1.2

    return R, T_comb, P_comb, T_SL, P_SL, g_0, m_dot_design, altitudes, gamma


def calculate_bell(P_design_bell, P_variation, P_comb, gamma, T_comb, R):
    ## Bell
    design_altitude_bell_ind = \
    np.where(abs(P_variation - P_design_bell) == np.min(abs(P_variation - P_design_bell)))[0][0]
    design_alt_bell = altitudes[design_altitude_bell_ind]

    A_ratio, Mj, Tj = area_ratio_calc(P_comb, P_design_bell, gamma, T_comb)
    v_e = exhaust_velocity(gamma, R, Tj, Mj)
    thrust_design = m_dot_design * v_e

    T_choked, rho_choked, c_choked, d_throat = choked_conditions(P_comb, T_comb, R, gamma, m_dot_design)

    A_ratio, Mj, Tj = area_ratio_calc(P_comb, P_design_bell, gamma, T_comb)
    A_exit = exit_area(A_ratio, d_throat)
    print(A_exit)
    bell_thrust = thrust_output(P_variation, P_comb, gamma, P_design_bell, thrust_design, A_exit)

    return bell_thrust, design_alt_bell, thrust_design


def calculate_spike(P_design_aero, P_variation, P_comb, gamma, T_comb, R):
    ## Aerospike

    design_altitude_aero_ind = \
    np.where(abs(P_variation - P_design_aero) == np.min(abs(P_variation - P_design_aero)))[0][0]
    design_alt_aero = altitudes[design_altitude_aero_ind]

    A_ratio, Mj, Tj = area_ratio_calc(P_comb, P_design_aero, gamma, T_comb)
    v_e = exhaust_velocity(gamma, R, Tj, Mj)
    thrust_design = m_dot_design * v_e
    T_choked, rho_choked, c_choked, d_throat = choked_conditions(P_comb, T_comb, R, gamma, m_dot_design)

    A_ratio, Mj, Tj = area_ratio_calc(P_comb, P_variation, gamma, T_comb)
    v_e = exhaust_velocity(gamma, R, Tj, Mj)

    aero_thrust = m_dot_design * v_e

    return aero_thrust, design_alt_aero, thrust_design


def pressure_altitude_map(altitudes, R):
    # read altitude density file
    # file = open('../auxiliary/Alt_Density_Temperature.txt')
    file = open('Alt_Density_Temperature.txt')
    reader = csv.reader(file, delimiter=' ', skipinitialspace=True)

    altitude = []
    density = []
    temperature = []
    for row in reader:
        # row = [alt,density,temperature]
        alt_column = float(row[0])
        density_column = float(row[1])
        temperature_column = float(row[2])

        altitude.append([alt_column])
        density.append([density_column])
        temperature.append([temperature_column])

    altitude_data = np.array([altitude])[0, :, 0] * 1e3
    density_data = np.array([density])[0, :, 0] * 1e3
    temperature_data = np.array([temperature])[0, :, 0]

    density = np.interp(altitudes, altitude_data, density_data)
    temperature = np.interp(altitudes, altitude_data, temperature_data)

    R = 287.058
    pressure = density * temperature * R

    ## For experimental case we require negative altitudes
    h = np.linspace(-25e3, 0, 625)
    P = 101325 * (1 - 2.25577 * 10 ** -5 * h) ** 5.25588

    altitudes = np.concatenate((h[:-1], altitudes))
    pressure = np.concatenate((P[:-1], pressure))

    return pressure, altitudes


def design_pressure(alt, P, altitudes):
    ind = np.where(np.abs(altitudes - alt) == np.min(np.abs(altitudes - alt)))[0][0]

    P_design = P[ind]

    return P_design


if __name__ == "__main__":

    condition = 'hot'
    P_amb = 101325
    experimental_design = True

    R, T_comb, P_comb, T_SL, P_SL, g_0, m_dot_design, altitudes, gamma = set_constants(condition)

    P_variation, altitudes = pressure_altitude_map(altitudes, R)

    if experimental_design == False:
        P_design_bell_alt = 14.6e3
        P_design_aero_alt = 14.6e3

        P_design_bell = design_pressure(P_design_bell_alt, P_variation, altitudes)
        bell_thrust, design_alt_bell, thrust_design_bell = calculate_bell(P_design_bell, P_variation, P_comb, gamma,
                                                                          T_comb, R)

        P_design_aero = design_pressure(P_design_aero_alt, P_variation, altitudes)
        aerospike_thrust, design_alt_aero, thrust_design_aero = calculate_spike(P_design_aero, P_variation, P_comb,
                                                                                gamma, T_comb, R)

        NPR_design_bell = P_comb / P_design_bell
        NPR_design_aero = P_comb / P_design_aero

        Isp_bell = bell_thrust / (m_dot_design * g_0)
        Isp_aero = aerospike_thrust / (m_dot_design * g_0)
        Isp_design_bell = thrust_design_bell / (m_dot_design * g_0)
        Isp_design_aero = thrust_design_aero / (m_dot_design * g_0)

        NPR_current = P_comb / P_variation

        NPR_ratio_bell = NPR_current / NPR_design_bell
        NPR_ratio_aero = NPR_current / NPR_design_aero

        Isp_scaling_bell = ((Isp_bell / Isp_design_bell) ** 2) ** (gamma / (gamma - 1))
        Isp_scaling_aero = ((Isp_aero / Isp_design_aero) ** 2) ** (gamma / (gamma - 1))

        plt.figure()
        plt.semilogx(NPR_current, Isp_bell, 'r')
        plt.semilogx(NPR_current, Isp_aero, 'k')
        plt.legend(['CD', 'Spike'])
        plt.ylabel('Isp (s)')
        plt.xlabel('NPR')

        plt.figure()
        plt.plot(altitudes, Isp_bell, 'r')
        plt.plot(altitudes, Isp_aero, 'k')
        plt.legend(['CD', 'Spike'])
        plt.ylabel('Isp (s)')
        plt.xlabel('Altitude (m)')

    else:
        NPR_design_test = 45
        throttle = np.array([0.1, 0.3, 0.5, 0.75, 1.0])
        # P_comb is combustion chamber pressure. Throttle is just a ratio
        NPR_throttles = (throttle * P_comb) / P_amb
        thrust_max = 3333.3
        Isp_max = 272.43

        # Randomly sets a nozzle pressure ratio of 45. And then just gets out a P_design?
        # Is P_design the pressure to set the exit area?
        P_design = P_comb / NPR_design_test

        #P_variation is the atmospheric pressure
        # Give me the indicy where the atmospheric - the 'design' is smallest. Where does the design pressure equal the atmospheric?
        P_ind = np.where(abs(P_variation - P_design) == np.min(abs(P_variation - P_design)))[0][0]

        condition = 'hot'

        P_design_bell_alt = altitudes[P_ind]
        P_design_aero_alt = altitudes[P_ind]

        # works out the design pressure of the bell?
        P_design_bell = design_pressure(P_design_bell_alt, P_variation, altitudes)
        bell_thrust, design_alt_bell, thrust_design_bell = calculate_bell(P_design_bell, P_variation, P_comb, gamma,
                                                                          T_comb, R)

        P_design_aero = design_pressure(P_design_aero_alt, P_variation, altitudes)
        aerospike_thrust, design_alt_aero, thrust_design_aero = calculate_spike(P_design_aero, P_variation, P_comb,
                                                                                gamma, T_comb, R)
        NPR_design_bell = P_comb / P_design_bell
        NPR_design_aero = P_comb / P_design_aero

        Isp_bell = bell_thrust / (m_dot_design * g_0)
        Isp_aero = aerospike_thrust / (m_dot_design * g_0)
        Isp_design_bell = thrust_design_bell / (m_dot_design * g_0)
        Isp_design_aero = thrust_design_aero / (m_dot_design * g_0)

        NPR_current = P_comb / P_variation

        m_dot_test = throttle * thrust_max / (Isp_max * g_0)
        bell_actual = np.zeros(5)
        aero_actual = np.zeros(5)

        Isp_bell_throttle = np.zeros(5)
        Isp_aero_throttle = np.zeros(5)
        for i in range(5):
            Isp_bell_throttle_ind = \
            np.where(abs(NPR_current - NPR_throttles[i]) == np.min(abs(NPR_current - NPR_throttles[i])))[0][0]
            Isp_bell_throttle[i] = Isp_bell[Isp_bell_throttle_ind]
            Isp_aero_throttle[i] = Isp_aero[Isp_bell_throttle_ind]

            bell_actual[i] = bell_thrust[Isp_bell_throttle_ind] / (m_dot_design) * m_dot_test[i]
            aero_actual[i] = aerospike_thrust[Isp_bell_throttle_ind] / (m_dot_design) * m_dot_test[i]

        difference_actual = (aero_actual - bell_actual) / bell_actual

        plt.figure()
        plt.plot(NPR_current, Isp_bell, 'r')
        plt.plot(NPR_current, Isp_aero, 'k')
        plt.legend(['CD', 'Spike'])
        plt.scatter(NPR_throttles, Isp_bell_throttle, c='r', marker='d')
        plt.scatter(NPR_throttles, Isp_aero_throttle, c='k', marker='d')
        plt.ylabel('Isp (s)')
        plt.xlabel('NPR')

        plt.show()
        plt.close()
        exit(0)

        plt.figure()
        plt.plot(altitudes, bell_thrust, 'r')
        plt.plot(altitudes, aerospike_thrust, 'k')
        plt.legend(['CD', 'Spike'])
        plt.ylabel('Thrust (N)')
        plt.xlabel('Altitude (m)')

        plt.show()
        plt.close()


        plt.figure()
        plt.plot(NPR_throttles, bell_actual, '.r')
        plt.plot(NPR_throttles, aero_actual, '.k')
        plt.legend(['CD', 'Spike'])
        plt.ylabel('Thrust (N)')
        plt.xlabel('NPR')

        plt.show()
        plt.close()