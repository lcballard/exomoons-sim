import numpy as np
from matplotlib import pyplot as plt
import time

# one dimensional model of planet and moon transiting sun

R_p = 1 * 6371  # earth radii to km
R_m = 0.2 * 6371  # earth radii to km
R_s = 1 * 696340 # solar radii to km
p_orbit_per = 1 # years DO NOT CHANGE UNIT (for using keplers 3rd)
p_orbit_e = 1   # eccentricity
m_orbit_per = 0.03 # years
m_orbit_e = 1
p_phase = 0
m_phase = 0
tmax = 1000 # number of days to simulate
n = tmax * 24 * 60 * 60 #timestep (in sec)
t = np.linspace(0,tmax,n)
pi = 3.1415


print("calculating orbits... ")
ti = time.time()
#calculate all planet positions
SMA_p = np.cbrt(p_orbit_per**2)
planet_x = SMA_p * np.sin((2*pi / (p_orbit_per*365)) * t + p_phase)
planet_x_dt = SMA_p * np.cos((2*pi / (p_orbit_per*365)) * t + p_phase)

#calculate all moon positions
SMA_m = np.cbrt(m_orbit_per**2)
moon_x = planet_x + SMA_m * np.sin((2*pi / (m_orbit_per*365)) * t + m_phase)
moon_x_dt = SMA_m * np.cos((2*pi / (m_orbit_per*365)) * t + m_phase)


#in_transitx = []
#in_transity = []
luminosity = np.ones(n) * 1 # 1 solar luminosity

#function for the area of intersection of 2 circles
def intersectionArea(X1, Y1, R1, X2, Y2, R2):
     
    # Calculate the euclidean distance
    # between the two points
    d = np.sqrt(((X2 - X1) * (X2 - X1)) + ((Y2 - Y1) * (Y2 - Y1)))
 
    if (d > R1 + R2):
        ans = 0
 
    elif (d <= (R1 - R2) and R1 >= R2):
        ans = pi * R2 * R2
 
    elif (d <= (R2 - R1) and R2 >= R1):
        ans = pi * R1 * R1
 
    else:
        alpha = np.arccos(((R1 * R1) + (d * d) - (R2 * R2)) / (2 * R1 * d)) * 2
        beta = np.arccos(((R2 * R2) + (d * d) - (R1 * R1)) / (2 * R2 * d)) * 2
         
        a1 = (0.5 * beta * R2 * R2 ) - (0.5 * R2 * R2 * np.sin(beta))
        a2 = (0.5 * alpha * R1 * R1) - (0.5 * R1 * R1 * np.sin(alpha))
        ans = a1 + a2
 
    return ans

print("... done!")
print(f"time: {((time.time() - ti)):.2f}s")

print("calculating transits... ")
ti = time.time()
for i in range(len(t)):
    tot_intersect = 0
    
    #check if planet is in front of sun at all
    if (planet_x[i]*149597870.691 + R_p > -R_s) and (planet_x[i]*149597870.691 - R_p < R_s):
        if planet_x_dt[i] < 0:
            #print(f"sun transits planet at t = {t[i]}")
            continue
        elif planet_x_dt[i] > 0:
            print(f"!! planet transits sun at t = {t[i]}")
            #in_transitx.append(t[i])
            #in_transity.append(planet_x[i])

            #print(r**2 * np.arccos((d**2+r**2-R**2)/(2*d*r)+0j) + R**2*np.arccos((d**2+R**2-r**2)/(2*d*R)+0j))

            planet_intersect = intersectionArea(planet_x[i]*149597870.691,0,R_p,0,0,R_s)
            #np.real(r**2 * np.arccos((d**2+r**2-R**2)/(2*d*r)+0j)) + np.real(R**2*np.arccos((d**2+R**2-r**2)/(2*d*R)+0j)) - 0.5*np.sqrt((-d+r+R)*(d+r-R)*(d-r+R)*(d+r+R))

            #print(area_intersect)
            #TODO: GET PERCENTAGE OF LUMINOSITY SUBTRACTION
            tot_intersect = planet_intersect

    #check if moon is in front of sun at all
    if (moon_x[i]*149597870.691 + R_m > -R_s) and (moon_x[i]*149597870.691 - R_m < R_s):
        print("moon transits sun!")
        area_intersect_m = intersectionArea(moon_x[i]*149597870.691,0,R_m,0,0,R_s)

        if (moon_x[i]*149597870.691 + R_m > -R_s) and (moon_x[i]*149597870.691 - R_m < R_s):
            moon_planet_overlap = intersectionArea(moon_x[i],0,R_m,planet_x[i],0,R_p)
            tot_intersect += area_intersect_m - moon_planet_overlap

    # set luminosity = % total luminosity uncovered
    if tot_intersect != 0:
        #print(f"{tot_intersect} = (moon: {area_intersect_m} - double: {moon_planet_overlap}) + planet: {planet_intersect}")
        #print(f"lumnosity factor: {1 - (tot_intersect)/(pi*R_s**2)}")
        luminosity[i] = luminosity[i] * (1 - (tot_intersect)/(pi*R_s**2))



    #DEBUG 
    if i % (0.1 * n) == 0:
        print(f"... finished step {i}/{n} [{i/n*100:.1f}%], at day {t[i]:.0f} ...")



#ax = plt.figure()
#plt.plot(t, planet_x)
#plt.xlabel("time(days)")
#plt.ylabel("position (AU)")
#plt.scatter(in_transitx, in_transity,c="k")
#plt.show()

print("... done!")
print(f"time: {((time.time() - ti)):.2f}s")
ax = plt.figure(figsize=(12,6))
plt.plot(t,luminosity)
plt.yticks([0.999999])
plt.xlabel("time (days)")
plt.ylabel("Luminosity ($M_\odot$)")
plt.show()


ax = plt.figure(figsize=(12,6))
plt.plot(moon_x,t)
plt.title("Moon position over time")
plt.xlabel("Distance from star (AU)")
plt.ylabel("time (days)")
plt.show()
