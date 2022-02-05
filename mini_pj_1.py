# Mini Project 1

import numpy as np
import matplotlib.pyplot as plt
import photovoltaic as pv


###################################################################################################################
# equation of time function
# input: julian day of the year
# output: equation of time
###################################################################################################################


def equation_of_time(day_no):
    B = 360.0 / 365.0 * (day_no - 81.0)
    EoT = 9.87 * pv.sind(2 * B) - 7.53 * pv.cosd(B) - 1.5 * pv.sind(B)
    # print('EoT', EoT)
    return EoT


###################################################################################################################
# time correction factor function
# input: equation of time, longitude, GMT offset of location
# output: time correction factor
###################################################################################################################
def time_correction(EoT, longitude, GMTOffset):
    LSTM = 15.0 * GMTOffset
    TimeCorrection = 4.0 * (longitude - LSTM) + EoT
    return TimeCorrection


###################################################################################################################
# local solar time function
# input: time correction, and local time
# output: local solar time
###################################################################################################################
def local_solar_time(TimeCorrection, H):
    local_solar_time = H + (TimeCorrection / 60)
    return local_solar_time


###################################################################################################################
# declination function
# input: julian day number
# output: declination of sun in degrees
###################################################################################################################
def declination(dayNo):
    x = (360 / 365) * (dayNo + 284)
    declination = 23.45 * pv.sind(x)
    return declination


###################################################################################################################
# elevation azimuth function
# input: declination of sun, latitude of location, and local solar time
# output: elevation and azimuth in degrees
###################################################################################################################
def elev_azi(declination, latitude, local_solar_time):
    hour_angle = 15.0 * (local_solar_time - 12.0)
    elevation = pv.arcsind(
        pv.sind(declination) * pv.sind(latitude) + pv.cosd(declination) * pv.cosd(latitude) * pv.cosd(hour_angle))
    azimuth = pv.arccosd(
        (pv.cosd(latitude) * pv.sind(declination) - pv.cosd(declination) * pv.sind(latitude) * pv.cosd(
            hour_angle)) / pv.cosd(elevation))
    # the multiplication by 1.0 causes a single value return for single inputs, otherwise it returns an array of one
    # element
    azimuth = np.where(hour_angle > 0, 360.0 - azimuth, azimuth) * 1.0
    return elevation, azimuth


###################################################################################################################
# sun position function
# input: julian day number, latitude, longitude, GMT offset of location, Hour, and Minute
# output: elevation and azimuth, position of the sun
###################################################################################################################
def sun_position(dayNo, latitude, longitude, GMTOffset, H, M):
    EoT = equation_of_time(dayNo)
    TimeCorrection = time_correction(EoT, longitude, GMTOffset)
    local_solar_time = H + (TimeCorrection + M) / 60.0
    elevation, azimuth = elev_azi(declination(dayNo), latitude, local_solar_time)
    return elevation, azimuth

#######################################################################################################################
# Problem 1
# read the given TMY data and plot the direct, diffuse, and global radiation for every hour of the year
# read CSV file and extract the needed columns of data
fname = '722970TYA.CSV'
station = np.genfromtxt(fname, max_rows=1, delimiter=",", usecols=(0))
location_name, location_state = np.genfromtxt(fname, max_rows=1, delimiter=",", usecols=(1, 2), dtype=str)
GHI, DNI, DHI = np.genfromtxt(fname, skip_header=2, delimiter=",", usecols=(4, 7, 10), unpack=True)
print('The data is from station number', station, 'at', location_name, '', location_state)

# plotting DNI for the year
plt.figure('DNI from TMY')
plt.title('DNI from TMY Long Beach, CA')
plt.xlabel('hours in the year')
plt.ylabel('direct normal irradiance')
plt.plot(DNI)

# plotting DHI for the year
plt.figure('DHI from TMY')
plt.title('DHI from TMY Long Beach, CA')
plt.xlabel('hours in the year')
plt.ylabel('direct normal irradiance')
plt.plot(DHI)

# plotting GHI for the year
plt.figure('GHI from TMY')
plt.title('GHI from TMY Long Beach, CA')
plt.xlabel('hours in the year')
plt.ylabel('global horizontal irradiance')
plt.plot(GHI)

######################################################################################################################
# Problem 2 and 2b
# calculate the position of the sun for every hour of the year, plot one day of this on polar plot
GMTOffset = -8  # Long Beach, CA is UTC-8 time zone
latitude = 33.8178  # for Long Beach Airport
longitude = -118.15167  # for Long Beach Airport
dayNo_plot = 120  # julian date of the year for the plot
M = 0  # position calculated on the hour
H = np.arange(0, 24, 1)  # hours of the day
days = np.arange(1, 366, 1)  # generate the days of the year
dayNo = np.repeat(days, 24)  # hours of the day, days of the year times
hours = np.arange(1, 25)  # Generate an array with the values 1 to 24
hour_of_day = np.tile(hours, 365)  # Generate an array with the hours 1 to 24 repeated 365 times

# calculate all sun positions every hour of the year
lb_declination = declination(dayNo)  # suns declination angle
lb_LST = local_solar_time(time_correction(equation_of_time(dayNo), longitude, GMTOffset), hour_of_day)
# 2 ###############################
elevation, azimuth = elev_azi(lb_declination, latitude, lb_LST)
###################################

# debug plot to check years worth of data
#theta_ = azimuth
#r_ = elevation
#plt.figure('Sun Path Over the Year')
#ax = plt.subplot(111, projection='polar')  # setup for sun path polar plot
#ax.plot(np.radians(theta_), 90.0 - r_, marker='.', linestyle='solid')  # setup for sun path polar plot
#ax.set_theta_zero_location('N')  # setup for sun path polar plot
#ax.set_theta_direction(-1)  # setup for sun path polar plot
#ax.set_rmax(90)  # setup for sun path polar plot
#ax.grid(True)  # setup for sun path polar plot
#ax.set_title("Path of the Sun at LB Daugherty Field, CA Over a Year", va='bottom')  # set title of plot
#plt.show()  # show the plot

# data to be used in julian day 120 plot
lb_declination_plot = declination(dayNo_plot)  # suns declination angle
lb_LST_plot = local_solar_time(time_correction(equation_of_time(dayNo_plot), longitude, GMTOffset), H)  # LST in Long Beach
elevation_plot, azimuth_plot = elev_azi(lb_declination_plot, latitude, lb_LST_plot)  # elevation, azimuth over course of the chosen day

# 2b ###########################
# setup for the polar plot
theta = azimuth_plot
r = elevation_plot
plt.figure('Sun Path')
ax = plt.subplot(111, projection='polar')  # setup for sun path polar plot
ax.plot(np.radians(theta), 90.0 - r, marker='.', linestyle='solid')  # setup for sun path polar plot
ax.set_theta_zero_location('N')  # setup for sun path polar plot
ax.set_theta_direction(-1)  # setup for sun path polar plot
ax.set_rmax(90)  # setup for sun path polar plot
ax.grid(True)  # setup for sun path polar plot
ax.set_title("Path of the Sun at LB Daugherty Field, CA on Day 120", va='bottom')  # set title of plot
##################################

#######################################################################################################################
# Problem 3
# calculate the normally incident solar radiation on a tilted surface for every hour of the year
beta = 33.8178  # tilt of the surface
phi = 180  # azimuth of the tilted surface
zenith = 90.0 - elevation  # angle of the sun from the vertical plane
#AM = 1/(pv.cosd(zenith))  # Air Mass calculation # debug statement results in nan values
AM = 1.5  # standard solar spectrum AM

# 3a #########################
S_incident = ((198 - (0.1 * beta))/180) * (1.353 * (0.7 ** (AM ** 0.678)))  # solar radiation incident on surface every hour
##############################

#print(S_incident)  # debug statement
# solar radiation normal to the surface every hour
S_module = S_incident * (pv.cosd(elevation) * pv.sind(beta) * pv.cosd(phi - azimuth) + pv.sind(elevation) * pv.cosd(beta))
I_G = S_module * 1.1  # total solar radiation on the surface, includes 10% extra for diffuse radiation every hour of the year
# loop to remove "negative power"
i = 0  # set up index
for x in I_G:
    if x < 0.0:
        I_G[i] = 0
        i += 1
    else:
        i += 1
i = 0  # set up index
I_G_year_total = 0  # placeholder for the for loop
for x in I_G:
    I_G_year_total += I_G[i]
    i += 1

# 3b ###########################
print('The total power density on the tilted surface for the year is', I_G_year_total, 'Wh/m^2')
################################

# 3c ###########################
plt.figure('PD hour of the year')
plt.title('Solar Radiation Power Density Yearly Total: 2684.86 Wh/m^2')
plt.xlabel('hours in the year')
plt.ylabel('Hourly Solar Radiation Power Density Wh/m^2')
plt.plot(I_G)
################################

#######################################################################################################################
# Problem 4
# Trace the path of the sun on a polar plot for the same time of day for every day of the year
# data to be used in same time plot
TOD = 16
days = np.arange(1, 366, 1)
lb_declination_hour = declination(days)  # suns declination angle
lb_LST_hour = local_solar_time(time_correction(equation_of_time(days), longitude, GMTOffset), TOD)  # LST in Long Beach
elevation_hour, azimuth_hour = elev_azi(lb_declination_hour, latitude, lb_LST_hour)  # elevation, azimuth over course of the chosen day
# setup for the polar plot
theta_hour = azimuth_hour
r_hour = elevation_hour
plt.figure('Sun Path at Same Time of Day')
ax = plt.subplot(111, projection='polar')  # setup for sun path polar plot
ax.plot(np.radians(theta_hour), 90.0 - r_hour, marker='.', linestyle='solid')  # setup for sun path polar plot
ax.set_theta_zero_location('N')  # setup for sun path polar plot
ax.set_theta_direction(-1)  # setup for sun path polar plot
ax.set_rmax(90)  # setup for sun path polar plot
ax.grid(True)  # setup for sun path polar plot
ax.set_title("Position of the Sun at LB Daugherty Field, CA on at 16:00 Each Day", va='bottom')  # set title of plot
plt.show()  # show the plot

