# Mini Project 3

import numpy as np
import matplotlib.pyplot as plt
import photovoltaic as pv

# mini project 2 values
J0 = 1.453958092629679e-13  # A/cm²
JL = 0.038909958438950246  # A/cm²
Area = 243.36  # cm²  # area of the solar cell
Rseries = 0.7522092113029638  # ohm.cm
voc = pv.cell.Voc(JL, J0)
ff_0 = 0.8027857996574461

# Read TMY data file, import by station name

fname = '722970TYA.CSV'
city = 'Long Beach, CA'
station, GMT_offset, latitude, longitude, altitude = np.genfromtxt(fname, max_rows=1, delimiter=",",
                                                                   usecols=(0, 3, 4, 5, 6))
ETR, GHI, DNI, DHI, ambient_temperature = np.genfromtxt(fname, skip_header=2, delimiter=",", usecols=(2, 4, 7, 10, 31),
                                                        unpack=True)

module_elevation = latitude
module_azimuth = 180  # point the module south
module_tilt = 0

print('station number: ', station)
print('City: ', city)
print('GMT: ', GMT_offset)
print('Latitude (degrees):', latitude)
print('Longitude (degrees):', longitude)
print('Altitude (m): ', altitude)

# details for the plots which makes string to label the plots
details4plot = 'lat: {:.1f}°, long: {:.1f}°, Module el: {}°, az: {}°'.format(
    latitude, longitude, module_elevation, module_azimuth)
print(details4plot)

# Part 1 ###############################################################################################################

# answer to question 1.a. ##############
I_cell = JL * 243.36  # (A)
# print(I_cell)  # debug statement
voltage = voc * 72  # 72 cells in series
# answer to question 1.b. ###############
P_module = I_cell * voltage * ff_0  # power output of a 72 cell module at the maximum power point
print('Power output of 72 cell module under AM 1.5 conditions (W):', P_module)

# Calculate the light intensity;
# Since JL from project 2 has been calculated at AM1.5G (1000W/m²),
# we need to divide by 1000 and multiply by the DHI and DNI to find JL for every hour of the year
# The unit of JL_module is A/cm²
JL_module = (DNI + DHI) * 0.001 * JL

# answer to question 2.a. #############
# This prints the current in mA/cm² from the module for every hour of the year
plt.figure('current_density')
plt.title('Current Density from the Module    AM1.5G Jsc = {:.2f} mA/cm²'.format(1e3 * JL))
plt.plot(1e3 * JL_module)
plt.xlabel('hour of the year')
plt.ylabel('Current Density in hour interval (mA/cm²)')
plt.show()

# We need to multiply by 1e4 to convert from 1/cm² to 1/m²
# and divide by 1e3 to convert from W to kW, resulting in a factor of 10
P_module = 10 * JL_module * voc * ff_0
# print('capacity factor = ', 100 * sum(P_module) / (10 * 8760 * JL * voc * ff_0), '%')
P_module = np.zeros(len(JL_module))  # Set up an array with zeros for the module power in each hour of the year
efficiency = np.zeros(len(JL_module))  # Set up an array with zeros for the module efficiency in each hour of the year
# Calculate the FF and efficiency for every hour the year
# We can used closed form equations for Rseries, but we need full IV curve for Part 2 anyway, so we will calcualte the IV curve (rather than just power) for every hour of the year
# idx is the looping index for calcualting the IV curve
for idx, JLm in enumerate(
        JL_module):  # Make a for loop; idx is indexing parameter for hours of the year, JLm is light generated current for the hour that we are calcualting
    Voc = pv.cell.Voc(JLm, J0)  # Calculate Voc for one hour in the year
    voltage = np.linspace(0, Voc,
                          300)  # Define 100 voltages to make an IV curve, going from V= 0 to V=v Voc for one hour of the year
    current = pv.cell.I_cell(voltage, JLm, J0)  # Calculate current for those 100 voltage points
    voltage_r = voltage - current * Rseries  # Calculate voltage impacted by Series Resistance for each of the 100 points in the IV curve
    Voc_r, Jsc_r, FF_r, Vmp_r, Jmp_r = pv.cell.cell_params(voltage_r,
                                                           current)  # Extract IV curve parameters for the one hour of the year
    P_module[
        idx] = Vmp_r * Jmp_r * 10  # Put the power from the module from that hour into an array called P_module. We ultiply by 1e4 to convert from W/cm² to W/m² and divide by 1e3 to convert to kW/m² from W/² and

P_yearly = sum(
    DNI + DHI) / 1e3  # Sum the total energy form the sun over the year. We divide by 1e3 to convert to kWh/m² from Wh/m²
P_yearly_module = sum(P_module)  # kWh/²  Sum of energy from the module over the year.
efficiency = P_yearly_module / P_yearly

# answers to question 2.b. ##############
print('\nUsing data from Long Beach, CA')
print('Energy from Sun over year {:.2f} kWh/m²'.format(sum(DNI + DHI) / 1000))
print('Energy from Module {:.2f} kWh/m²'.format(sum(P_module)))
print('Efficiency of the Module (%):', efficiency)
# plot the power density output
plt.figure('power_density')
plt.title('Power Density from the Module    Total Energy = {:.2f} kWh/m²'.format(P_yearly_module))
plt.plot(P_module)
plt.xlabel('hour of the year')
plt.ylabel('Power Density in hour interval (kWh/m²)')
plt.show()

# Part 2 ###############################################################################################################
# answers to question 1 are the plots output in first loop and question 2 is second loop #######
shading = 0.50
number_series = 8  # Number UNSHADED in series
shading_conc = (number_series + shading) / (number_series + 1)
number_points = 100  # Number of points in each IV curve
Vbr = [-5, -0.5]
for x in Vbr:
    # First we will calculate individual IV curves

    current = np.linspace(0, JL, number_points)

    # IV curve for unshaded solar cell
    voltage = 0.0248 * np.log((
                                      JL - 0.99999999999 * current) / J0)  # Put the .99999 in so we don't try to take a log of zero (we calcualte Voc seperately)
    voltage_r = voltage - current * Rseries
    Voc, Jsc, FF, Vmp, Jmp = pv.cell.cell_params(voltage_r, current)  # Extract IV curve parameters
    # Plot the  IV curves
    plt.plot(voltage, current)
    plt.plot(voltage_r, current)
    plt.scatter(Vmp, Jmp, label='Maximum Power Point')
    plt.text(Vmp, Jmp, ' Maximum \nPower Point')
    plt.ylim(0)
    plt.twinx()
    plt.plot(voltage_r, voltage_r * current, c='red')
    plt.ylim(0)
    plt.show()

    # IV curve for shaded solar cell
    JL_shaded = JL * shading
    current_no_zero = np.maximum((JL_shaded - current), 1e-30)
    voltage_positive_shaded = 0.0248 * np.log((current_no_zero / J0) + 1)
    voltage_r_positive = voltage_positive_shaded - current_no_zero * Rseries
    reverse_bias_check = np.where(JL_shaded - current > 0, 0, 1)
    voltage_shaded = (1 - reverse_bias_check) * voltage_r_positive + reverse_bias_check * x

    # Combined IV curve
    voltage_total = number_series * voltage_r + voltage_shaded
    Voc_total, Jsc_total, FF_total, Vmp_total, Jmp_total = pv.cell.cell_params(voltage_total,
                                                                               current)  # Extract IV curve parameters

    # Plot the Unshaded singe cell, shaded, and total  IV curves
    plt.plot(voltage_r, current)
    plt.plot(voltage_shaded, current)
    plt.plot(voltage_total, current)
    plt.xlabel('Voltage (V)')
    plt.ylabel('Current Density (A/cm²)')
    plt.ylim(0)
    plt.show()

    # Plot the total IV curve and Maximum Power point
    label = 'Project 3 Part 2 Module Shaded Series IV Curves '
    plt.title('IV Curves for Shaded Series Cells with {:.0f} in series'.format(number_series + 1))
    # plt.figure(label)
    plt.plot(voltage_total, current)
    plt.scatter(Vmp_total, Jmp_total, label='Maximum Power Point')
    plt.text(Vmp_total, Jmp_total, ' Maximum \nPower Point')
    plt.xlim(-1)
    plt.xlabel('Voltage (V)')
    plt.ylabel('Current Density (A/cm²)')
    # plt.legend()
    plt.ylim(0)
    plt.twinx()
    plt.plot(voltage_total, voltage_total * current, c='red')
    plt.ylim(0)
    plt.ylabel('Power Density(W/cm²)')
    plt.show()

P_module_unshaded = np.zeros(len(JL_module))  # Set up an array with zeros for the module power in each hour of the year
P_module_shaded = np.zeros(len(JL_module))  # Set up an array with zeros for the module power in each hour of the year
for idx, JLm in enumerate(
        JL_module):  # Make a for loop; idx is indexing parameter for hours of the year, JLm is light generated current for the hour that we are calcualting
    if JLm > 0:
        current = np.linspace(0, JLm, number_points)
        # IV curve for unshaded solar cell
        voltage = 0.0248 * np.log((
                                              JLm - 0.99999999999 * current) / J0)  # Put the .99999 in so we don't try to take a log of zero (we calcualte Voc seperately)
        voltage_r = voltage - current * Rseries
        Voc, Jsc, FF, Vmp, Jmp = pv.cell.cell_params(voltage_r, current)  # Extract IV curve parameters
        P_module_unshaded[
            idx] = Vmp * Jmp * 10  # Put the power from the module from that hour into an array called P_module. We ultiply by 1e4 to convert from W/cm² to W/m² and divide by 1e3 to convert to kW/m² from W/² and
        # IV curve for shaded solar cell
        JL_shaded = JLm * shading
        current_no_zero = np.maximum((JL_shaded - current), 1e-30)
        voltage_positive_shaded = 0.0248 * np.log((current_no_zero / J0) + 1)
        voltage_r_positive = voltage_positive_shaded - current_no_zero * Rseries
        reverse_bias_check = np.where(JL_shaded - current > 0, 0, 1)
        voltage_shaded = (1 - reverse_bias_check) * voltage_r_positive + reverse_bias_check * -0.5
        # Combined IV curve
        voltage_total = number_series * voltage_r + voltage_shaded
        Voc_total, Jsc_total, FF_total, Vmp_total, Jmp_total = pv.cell.cell_params(voltage_total,
                                                                                   current)  # Extract IV curve parameters
        P_module_shaded[idx] = (Vmp_total * Jmp_total * 10) / (
                    number_series + 1)  # Put the power from the module from that hour into an array called P_module. We ultiply by 1e4 to convert from W/cm² to W/m² and divide by 1e3 to convert to kW/m² from W/² and
    else:
        P_module_unshaded[idx] = 0
        P_module_shaded[idx] = 0

P_yearly_module_unshaded = sum(P_module_unshaded)  # kWh/²  Sum of energy from the module over the year.
P_yearly_module_shaded = sum(P_module_shaded)  # kWh/²  Sum of energy from the module over the year.
efficiency_unshaded = P_yearly_module_unshaded / P_yearly
efficiency_shaded = P_yearly_module_shaded / P_yearly
print('Energy from Sun over year {:.2f} kWh/m²'.format(sum(DNI + DHI) / 1000))
# print('Energy from Unshaded Module {:.2f} kWh/m²'.format(sum(P_module_unshaded)))
# print('Efficiency of the Unshaded Module (%):', efficiency_unshaded)
print('Energy from Shaded Module with bypass diode {:.2f} kWh/m²'.format(sum(P_module_shaded)))
print('Efficiency of the Shaded Module (%):', efficiency_shaded)

# question 3 ########################################
P_module_part_shaded = np.zeros(
    len(JL_module))  # Set up an array with zeros for the module power in each hour of the year
for idx, JLm in enumerate(
        JL_module):  # Make a for loop; idx is indexing parameter for hours of the year, JLm is light generated current for the hour that we are calcualting
    if JLm > 0:
        current = np.linspace(0, JLm, number_points)
        # IV curve for unshaded solar cell
        voltage = 0.0248 * np.log((
                                              JLm - 0.99999999999 * current) / J0)  # Put the .99999 in so we don't try to take a log of zero (we calcualte Voc seperately)
        voltage_r = voltage - current * Rseries
        Voc, Jsc, FF, Vmp, Jmp = pv.cell.cell_params(voltage_r, current)  # Extract IV curve parameters
        P_module_unshaded[
            idx] = Vmp * Jmp * 10  # Put the power from the module from that hour into an array called P_module. We ultiply by 1e4 to convert from W/cm² to W/m² and divide by 1e3 to convert to kW/m² from W/² and
        # IV curve for shaded solar cell
        if idx % 24 in (7, 8, 9):
            JL_shaded = JLm * shading
        else:
            JL_shaded = JLm
        current_no_zero = np.maximum((JL_shaded - current), 1e-30)
        voltage_positive_shaded = 0.0248 * np.log((current_no_zero / J0) + 1)
        voltage_r_positive = voltage_positive_shaded - current_no_zero * Rseries
        reverse_bias_check = np.where(JL_shaded - current > 0, 0, 1)
        voltage_shaded = (1 - reverse_bias_check) * voltage_r_positive + reverse_bias_check * -0.5
        # Combined IV curve
        voltage_total = number_series * voltage_r + voltage_shaded
        Voc_total, Jsc_total, FF_total, Vmp_total, Jmp_total = pv.cell.cell_params(voltage_total,
                                                                                   current)  # Extract IV curve parameters
        P_module_part_shaded[idx] = (Vmp_total * Jmp_total * 10) / (
                    number_series + 1)  # Put the power from the module from that hour into an array called P_module. We ultiply by 1e4 to convert from W/cm² to W/m² and divide by 1e3 to convert to kW/m² from W/² and
    else:
        P_module_unshaded[idx] = 0
        P_module_part_shaded[idx] = 0
P_yearly_module_part_shaded = sum(P_module_part_shaded)  # kWh/m²  Sum of energy from the module over the year.
efficiency_part_shaded = P_yearly_module_part_shaded / P_yearly
print('Energy from Partly Shaded Module with bypass diode {:.2f} kWh/m²'.format(sum(P_module_part_shaded)))
print('Efficiency of the Shaded Module (%):', efficiency_part_shaded)

# question 4 #############################
print('Energy from Module at Long Beach, CA {:.2f} kWh/m²'.format(sum(P_module)))
print('Energy from Partly Shaded Module with bypass diode at Long Beach, CA {:.2f} kWh/m²'.format(sum(P_module_part_shaded)))
print('Ratio of Power due to Partly Shaded Panel {:.2f} kWh/m²'.format(sum(P_module_part_shaded) / sum(P_module)))
