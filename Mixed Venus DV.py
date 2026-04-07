import pykep as pk
from pykep import DAY2SEC
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# Anti-clockwise --------------------------------------------------------------
dti_depart_ACW = pd.date_range("2026-03-01", periods=900, freq="D")
dti_arrive_ACW = pd.date_range("2027-07-01", periods=600, freq="D")

earth = pk.planet.jpl_lp('earth')
venus = pk.planet.jpl_lp('venus')

DeltaVACW = np.full((len(dti_depart_ACW), len(dti_arrive_ACW)), np.nan)

for i, dep_date in enumerate(dti_depart_ACW):
    td = pk.epoch_from_string(f'{dep_date},00:00:00')
    rE, vE = earth.eph(td)
    rE = np.array(rE)
    vE = np.array(vE)
    for j, arr_date in enumerate(dti_arrive_ACW):
        ta = pk.epoch_from_string(f'{arr_date},00:00:00')
        dt = (ta.mjd2000 - td.mjd2000) * DAY2SEC
        rV, vV = venus.eph(ta)
        rV = np.array(rV)
        vV = np.array(vV)
        if 1 * DAY2SEC <= dt <= 500 * DAY2SEC:
            l = pk.lambert_problem(r1=rE, r2=rV, tof=dt, mu=pk.MU_SUN, cw = False, max_revs=0)
            v1 = np.array(l.get_v1()[0])
            v2 = np.array(l.get_v2()[0])
            dv_depart = np.linalg.norm(v1 - vE) / 1000.0
            dv_arrive = np.linalg.norm(v2 - vV) / 1000.0
            DeltaVACW[i, j] = dv_depart + dv_arrive

# Clockwise -------------------------------------------------------------------
dti_depart_CW = pd.date_range("2027-01-01", periods=500, freq="D")
dti_arrive_CW = pd.date_range("2027-08-01", periods=300, freq="D")

DeltaVCW = np.full((len(dti_depart_CW), len(dti_arrive_CW)), np.nan)


for i, dep_date in enumerate(dti_depart_CW):
    td = pk.epoch_from_string(f'{dep_date},00:00:00')
    rE, vE = earth.eph(td)
    rE = np.array(rE)
    vE = np.array(vE)
    for j, arr_date in enumerate(dti_arrive_CW):
        ta = pk.epoch_from_string(f'{arr_date},00:00:00')
        dt = (ta.mjd2000 - td.mjd2000) * DAY2SEC
        rV, vV = venus.eph(ta)
        rV = np.array(rV)
        vV = np.array(vV)
        if 1 * DAY2SEC <= dt <= 500 * DAY2SEC:
            l = pk.lambert_problem(r1=rE, r2=rV, tof=dt, mu=pk.MU_SUN, cw = True, max_revs=0)
            v1 = np.array(l.get_v1()[0])
            v2 = np.array(l.get_v2()[0])
            dv_depart = np.linalg.norm(v1 - vE) / 1000.0
            dv_arrive = np.linalg.norm(v2 - vV) / 1000.0
            DeltaVCW[i, j] = dv_depart + dv_arrive
            
DeltaVACW[DeltaVACW > 20] = np.nan
DeltaVCW[DeltaVCW > 90] = np.nan

fig, ax = plt.subplots(figsize=(12, 8))

# Meshgrid for anti-clockwise
X_ACW, Y_ACW = np.meshgrid(mdates.date2num(dti_depart_ACW), 
                           mdates.date2num(dti_arrive_ACW), indexing='ij')
# Meshgrid for clockwise
X_CW, Y_CW = np.meshgrid(mdates.date2num(dti_depart_CW), 
                         mdates.date2num(dti_arrive_CW), indexing='ij')

# Plots for both
csACW = ax.contourf(X_ACW, Y_ACW, DeltaVACW, levels=100, cmap='jet') 
#, linewidths=0.8 ## If ax.contour( is used then this can be included
ax.text(0.3, 0.55, 'Anti-clockwise', transform=ax.transAxes, fontsize=14)

csCW = ax.contourf(X_CW, Y_CW, DeltaVCW, levels=100, cmap='jet', alpha=0.6)
ax.text(0.65, 0.25, 'Clockwise', transform=ax.transAxes, fontsize=14)


ax.set_xlabel("Departure Date")
ax.set_ylabel("Arrival Date")
ax.xaxis_date()
ax.yaxis_date()
fig.autofmt_xdate()

# Colourbars
cb1 = fig.colorbar(csACW, ax=ax, pad=0.02, label='Anti-clockwise Total $\Delta v$ (km/s)')
for line in cb1.lines: 
   line.set_linewidth(10) # Gives the colourbar the rainbow effect
cb2 = fig.colorbar(csCW, ax=ax, pad=0.1, label='Clockwise Total $\Delta v$ (km/s)')
for line in cb2.lines: 
   line.set_linewidth(10) #Gives the colourbar the rainbow effect

plt.grid(alpha=0.3)
plt.title(r'$\Delta v$ requirements in km/s for orbital transfers from Earth to Venus')
plt.show()
