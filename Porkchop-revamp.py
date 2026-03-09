# Imports
import pykep as pk
from pykep import DAY2SEC
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

#Generate date ranges for departure and arrival 
dti_depart = pd.date_range("2026-03-01", periods=900, freq="D")
dti_arrive = pd.date_range("2027-07-01", periods=600, freq="D")

#Ascertain information for Earth and Venus 
earth = pk.planet.jpl_lp('earth')
venus = pk.planet.jpl_lp('venus')

AU_e = 1.496e+11 #m
AU_v = 1.08159e+11 #m
#Create iterated list of DeltaV for each date
DeltaV = np.full((len(dti_depart), len(dti_arrive)), np.nan)

mu_sun = 1.32712440018e+20
mu_earth = 398600441800000
mu_venus = 324859000000000

orb_rad_e = 7015800.0000000009
orb_rad_v = 6657200.0000000009

orbital_vE = np.sqrt(mu_sun/AU_e)
orbital_vV = np.sqrt(mu_sun/AU_v)

#dV to leave Earth
vb1 = np.sqrt(2 * mu_sun * ((1/AU_e) - (1/(AU_e + AU_v))))
h_tran_earth = (vb1 - orbital_vE) / 1000 #km/s

vb2 = np.sqrt(2 * mu_sun * ((1/AU_v) - (1/(AU_e + AU_v))))

vperiapsis_v = np.sqrt(((vb2 - orbital_vV) ** 2) + ((2 * mu_venus) / orb_rad_v))
vorbit_v = np.sqrt((mu_venus / orb_rad_v))
vcapture_v = (vperiapsis_v - vorbit_v) / 1000 #km/s
hohmann_dv = vcapture_v + h_tran_earth


for i, dep_date in enumerate(dti_depart):
    td = pk.epoch_from_string(f'{dep_date},00:00:00')
    rE, vE = earth.eph(td)
    rE = np.array(rE)
    vE = np.array(vE)
    for j, arr_date in enumerate(dti_arrive):
        ta = pk.epoch_from_string(f'{arr_date},00:00:00')
        dt = (ta.mjd2000 - td.mjd2000) * DAY2SEC
        rV, vV = venus.eph(ta)
        rV = np.array(rV)
        vV = np.array(vV)
        if 40 * DAY2SEC <= dt <= 500 * DAY2SEC:
            l = pk.lambert_problem(r1=rE, r2=rV, tof=dt, mu=pk.MU_SUN, cw = False, max_revs=0)
            v1 = np.array(l.get_v1()[0])
            v2 = np.array(l.get_v2()[0])
            dv_depart = np.linalg.norm(v1 - vE) / 1000.0
            dv_arrive = np.linalg.norm(v2 - vV) / 1000.0
            DeltaV[i, j] = dv_depart + dv_arrive
            
print(np.nanmin(DeltaV))
print(dv_depart)
print(dv_arrive)

#Create and mask an array of invalid DeltaV values to limit the porkchop plot
DeltaV[DeltaV > 20] = np.nan
masked_DeltaV = np.ma.masked_invalid(DeltaV)
#Convert epoch to dates
dep_num = mdates.date2num(dti_depart.to_pydatetime())
arr_num = mdates.date2num(dti_arrive.to_pydatetime())

# Use indexing='ij' to match matrix orientation
X, Y = np.meshgrid(dep_num, arr_num, indexing='ij')

fig, ax = plt.subplots(figsize=(10, 8))

cs = ax.contour(X, Y, DeltaV, levels=100, cmap='jet', linewidths=0.3)

ax.set_xlabel("Departure Date")
ax.set_ylabel("Arrival Date")
#Tell python to read numbers as dates
ax.xaxis_date()
ax.yaxis_date()

cb = fig.colorbar(cs, ax=ax, orientation='vertical', label=r"Total $\Delta v$ (km/s)")
for line in cb.lines: 
   line.set_linewidth(10) #Gives the colourbar the rainbow effect

fig.autofmt_xdate() # Autoformats the dates
plt.grid(alpha=0.2)
plt.title(r'Porkchop plot measuring $\Delta v$ in km/s for transfers from Earth to Venus')
plt.show()