# this program finds the declination angle and elevation of the sun at solar noon at specific locations and days of
# the year
import numpy as np

# constants
# julian dates
# January 29th, April 21st, August 23rd, September 21st
days = [29, 111, 235, 264]
# latitudes
# Phoenix, AZ  San Juan, PR  Nome, AK  New York, NY
lats = [33.4485, 18.4663, 64.5023, 40.7167]

# empty placeholders
elevation_sn = []  # list for elevations
declination = []  # list for declinations

# declination for specified julian dates
for i in range(len(days)):
    x = (360 / 365) * (days[i] + 284)
    delta = 23.45 * (np.sin(np.radians(x)))
    declination.append(delta)
print("Declination list:", declination)

# sun elevation for specified locations on specified julian dates at solar noon
for j in range(len(lats)):
    alpha = 90 - lats[j] + declination[j]
    elevation_sn.append(alpha)
print("Elevation list:", elevation_sn)


