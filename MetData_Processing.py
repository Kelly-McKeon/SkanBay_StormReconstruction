#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 17 16:57:26 2025

@author: kellymckeon
"""

import numpy as np
import os
import matplotlib.pyplot as plt
import datetime
import pandas as pd
from scipy import stats

mac = r"/Users/kellymckeon/Library/CloudStorage/OneDrive-WoodsHoleOceanographicInstitution/Documents/WHOI/Alaska/"
pc = r"\Users\kelly\OneDrive - Woods Hole Oceanographic Institution\Documents\WHOI\Alaska"

# CHANGE DIRECTORY TO WHICHEVER COMPUTER YOU ARE ON
os.chdir(mac)

#%% Load Instrument Data

metData = pd.read_csv(r"Codes/Instruments/MetData.csv", index_col="Time", parse_dates=True)
tiltData = pd.read_csv(r"Codes/Instruments/Tiltmeter.csv", index_col='ISO 8601 Time', parse_dates=True)
hoboData = pd.read_csv(r"Codes/Instruments/HOBO.csv",index_col="Time", parse_dates=True)
waveData_in = pd.read_csv(r"Codes/Instruments/Inner_Buoy.csv")
waveData_out = pd.read_csv(r"Codes/Instruments/Outer_Buoy.csv")

#wave data in epoch time so convert it to datetime
waveData_in['datetime'] = pd.to_datetime(waveData_in['time'], unit='s')
waveData_in.set_index('datetime')

waveData_out['datetime'] = pd.to_datetime(waveData_out['time'], unit='s')
waveData_out.set_index('datetime')

#get max windspeed and mean wind direction and mean temperature
metMax = metData.resample('D').max()
metMean = metData.resample('D').mean()
metMin = metData.resample('D').min()
maxDir_ind = metData.resample('D')[['Wspd_kts']].idxmax()['Wspd_kts'].dropna()
metMax['Wdir_degT'] = metData['Wdir_degT'].loc[maxDir_ind].to_numpy()

#get max current and mean direction
tilt =  tiltData.resample('H').mean()
tiltMax = tilt.resample('D').max()
tiltMean = tiltData.resample('D').mean()

#get max wave height and mean 
waveMax_in = waveData_in.resample('D', on='datetime').max()
waveMean_in = waveData_in.resample('D', on='datetime').mean()

waveMax_out = waveData_out.resample('D', on='datetime').max()
waveMean_out = waveData_out.resample('D', on='datetime').mean()


#get max water depth
hoboMax = hoboData.resample('D').max()
hoboMean = hoboData.resample('D').mean()


#%% Plotting Maximums
figs,axs = plt.subplots(4,1,figsize=(12,20),dpi=150)
ax0=axs[0]
ax1=axs[1]
ax2=axs[2]
ax3=axs[3]

#plot winds

#speeds
ax00 = ax0.twinx()
p1=ax0.plot(metMax["Wspd_kts"], color='blue', label = 'speed')
ax0.set_ylabel('Max Daily WindSpeed (kts)')

#direction
p2=ax00.plot(metMean['Wdir_degT'], color = 'orange', label = 'direction')
ax00.set_ylabel('Mean Wind Direction (deg T)')
ax00.set(ylim=(0,360))
leg = p1 + p2
lab = [l.get_label() for l in leg]
ax0.legend(leg,lab,loc=0)
ax0.set_title('Winds')



#plot buoys

#pressure on buoy plot
ax01=ax1.twinx()
p3=ax01.plot(metMin["BP_mbar"], label='pressure', color='black')
ax01.set_ylabel('Min Daily Atmospheric Pressure (mbar)')
ax01.set(ylim=(900,1050))


#buoy heights
p4=ax1.plot(waveMax_in["height"],label='inner buoy')
p5=ax1.plot(waveMax_out["height"], label = 'outer buoy')
ax1.set_ylabel('Max Daily Wave Height (m)')
leg = p3 + p4 + p5
lab = [l.get_label() for l in leg]
ax1.legend(leg,lab,loc=0)
ax1.set_title('Buoys')

#plot WL
ax2.plot(hoboMax['WL']/3.281) #convert ft to m
ax2.set_ylabel('Max Daily Water Level at Inlet (m)')
ax2.set_title('Water Level')

#plot tiltmeter
ax03=ax3.twinx()

#speed
p6=ax3.plot(tiltMax['Speed (cm/s)'], color = 'blue', label = 'speed')
ax3.set_ylabel('Current Speed(cm/s')

p7=ax03.plot(tiltMean['Heading (degrees)'], color = 'orange', label = 'heading')
ax03.set_ylabel('Current Heading (degrees)')
leg = p6 + p7
lab = [l.get_label() for l in leg]
ax3.legend(leg,lab,loc=0)
ax3.set_title('Current')

#%% Plotting Means
del figs,axs

figs,axs = plt.subplots(4,1,figsize=(12,20),dpi=150)
ax0=axs[0]
ax1=axs[1]
ax2=axs[2]
ax3=axs[3]

#plot winds

#speeds
ax00 = ax0.twinx()
p1=ax0.plot(metMean["Wspd_kts"], color='blue', label = 'speed')
ax0.set_ylabel('Mean Daily WindSpeed (kts)')

#direction
p2=ax00.plot(metMean['Wdir_degT'], color = 'orange', label = 'direction')
ax00.set_ylabel('Mean Wind Direction (deg T)')
ax00.set(ylim=(0,360))
leg = p1 + p2
lab = [l.get_label() for l in leg]
ax0.legend(leg,lab,loc=0)
ax0.set_title('Winds')



#plot buoys

#pressure on buoy plot
ax01=ax1.twinx()
p3=ax01.plot(metMean["BP_mbar"], label='pressure', color='black')
ax01.set_ylabel('Mean Daily Atmospheric Pressure (mbar)')
ax01.set(ylim=(900,1050))


#buoy heights
p4=ax1.plot(waveMean_in["height"],label='inner buoy')
p5=ax1.plot(waveMean_out["height"], label = 'outer buoy')
ax1.set_ylabel('Mean Daily Wave Height (m)')
leg = p3 + p4 + p5
lab = [l.get_label() for l in leg]
ax1.legend(leg,lab,loc=0)
ax1.set_title('Buoys')

#plot WL
ax2.plot(hoboMean['WL']/3.281) #convert ft to m
ax2.set_ylabel('Mean Daily Water Level at Inlet (m)')
ax2.set_title('Water Level')

#plot tiltmeter
ax03=ax3.twinx()

#speed
p6=ax3.plot(tiltMean['Speed (cm/s)'], color = 'blue', label = 'speed')
ax3.set_ylabel('Current Speed(cm/s')

p7=ax03.plot(tiltMean['Heading (degrees)'], color = 'orange', label = 'heading')
ax03.set_ylabel('Current Heading (degrees)')
leg = p6 + p7
lab = [l.get_label() for l in leg]
ax3.legend(leg,lab,loc=0)
ax3.set_title('Current')

#%% Plot Just Maximum Wind and Waves During October Event

figs,axs = plt.subplots(4,1,dpi=150, constrained_layout=True)
ax0=axs[0]
ax1=axs[1]
ax2=axs[2]
ax3=axs[3]

#plot winds

#speeds
ax00 = ax0.twinx()
p1=ax0.plot(metData["Wspd_kts"], color='tab:orange', label = 'speed')
ax0.set_ylabel('Wind Speed (kts)')

#direction
p2=ax00.plot(metData['Wdir_degT'], color = 'tab:blue', label = 'direction')
ax00.axhspan(300,360, color='k', alpha=0.1, label='300 degT')
ax00.set_ylabel('Wind Direction (deg T)')
ax00.set(ylim=(0,360))
leg = p1 + p2
lab = [l.get_label() for l in leg]
ax0.legend(leg,lab,loc=0)
ax0.grid(color='lightgray', linestyle=':', linewidth=0.7, axis='x', which='major')
ax0.set_xlim(datetime.date(2021,10,18), datetime.date(2021,10,24))
# ax0.set_title('Winds')



#plot buoys

#pressure on buoy plot
ax01=ax1.twinx()
p3=ax01.plot(metData["BP_mbar"], label='pressure', color='black')
ax01.set_ylabel('Atmospheric Pressure (mbar)')
ax01.set(ylim=(900,1050))
ax1.grid(color='lightgray', linestyle=':', linewidth=0.7, axis='x', which='major')
ax1.set_xlim(datetime.date(2021,10,18), datetime.date(2021,10,24))


#buoy heights
p4=ax1.plot(waveData_in['datetime'], waveData_in["height"],label='inner buoy')
p5=ax1.plot(waveData_out['datetime'], waveData_out["height"], label = 'outer buoy')
ax1.set_ylabel('Wave Height (m)')
leg = p3 + p4 + p5
lab = [l.get_label() for l in leg]
ax1.legend(leg,lab,loc=0)
ax1.set_xlim(datetime.date(2021,10,18), datetime.date(2021,10,24))


#plot WL
ax2.plot(hoboData['WL']/3.281) #convert ft to m
ax2.set_ylabel('Water Level (m)')
ax2.grid(color='lightgray', linestyle=':', linewidth=0.7, axis='both', which='major')
ax2.set_xlim(datetime.date(2021,10,18), datetime.date(2021,10,24))
# ax2.set_title('Water Level')

#plot tiltmeter
ax03=ax3.twinx()

#speed
p6=ax3.plot(tilt['Speed (cm/s)'], color = 'tab:blue', label = 'speed')
ax3.set_ylabel('Current Speed (cm/s)')
ax3.grid(color='lightgray', linestyle=':', linewidth=0.7, axis='x', which='major')
ax3.set_xlim(datetime.date(2021,10,19), datetime.date(2021,10,24))
ax3.set_ylim(30,50)

p7=ax03.plot(tilt['Heading (degrees)'], color = 'tab:orange', label = 'heading')
ax03.set_ylabel('Current Heading (degrees)')
ax03.set_ylim(200,360)
leg = p6 + p7
lab = [l.get_label() for l in leg]
ax3.legend(leg,lab,loc=0)
ax3.set_xlim(datetime.date(2021,10,18), datetime.date(2021,10,24))


#%% Load ASOS Airport data 

## downloaded from https://mesonet.agron.iastate.edu/request/download.phtml?network=AK_ASOS
## hourly wind speed and direction data since 1946

#timestamp is column 0
#wind direction (deg true north) is column 4 ('drct')
#wind speed (knots) is column 5 ('sknt')
#slp is column 6 ('mslp')
#peak wind gust is column 9 ('peak_wind_gust')
#peak wind gust direction is column 10 ('peak_wind_direction')
#peak wind time is column 11 ('peak_wind_time')


airData = pd.read_csv(r"Codes/Instruments/ASOS_Airport_Data.csv", index_col='time',parse_dates=True, date_format='%m/%d/%Y %H:%M')
# airData['datetime'] = pd.to_datetime(airData['time'], format='%m/%d/%Y %H:%M')

#get max and mean wind speeds and associated directions
speed_ind = airData.resample('D')[['sknt']].idxmax()['sknt'].dropna() #get indexes of daily maximum speeds
airSpeedMax = airData['sknt'].loc[speed_ind].to_numpy() #get daily max wind speed
airSpeedMean = airData['sknt'].resample('D').mean().to_numpy() #get daily average wind speed
airDirMax = airData['drct'].loc[speed_ind].to_numpy() #get wind direction during max speed
airDirMean = airData['drct'].resample('D').mean().to_numpy() #get average daily wind direction
windMax_time = speed_ind.index.to_series() #convert windspeed max times to plottable series
windMean_time = airData.index.to_series() #convert windspeed mean times to plottable series
windMean_time = windMean_time.resample('D').min()



#%% Load NDBC Wave Data

# NDBC Buoy 46073 205NM WNW of Dutch Harbor 
# data downloaded from https://portal.aoos.org/#metadata/41997/station/data
#times are in ISO8601 UTC
#time is column 0, height is 1, quality control is 2
#units are m, deg True 


#heights
ndbc_waveHgt = pd.read_csv('Codes/Instruments/NDBC_WaveHgt.csv', index_col='time', parse_dates=True)
height_ind = ndbc_waveHgt.resample('D')[['height']].idxmax()['height'].dropna()
ndbcHeightMax = ndbc_waveHgt['height'].loc[height_ind].to_numpy()
ndbcHeightMean = ndbc_waveHgt['height'].resample('D').mean().to_numpy()

#directions
ndbc_waveDir = pd.read_csv('Codes/Instruments/NDBC_WaveDir.csv', index_col='time', parse_dates=True)
ndbcDirMax = ndbc_waveDir['dir'].loc[height_ind].to_numpy()
ndbcDirMean = ndbc_waveDir['dir'].resample('D').mean().to_numpy()

ndbc_time = height_ind.index.to_series()



#%% Compare Skan Bay Met Data to Dutch Harbor Airport Data

fig, ax = plt.subplots(4,1, dpi=150, constrained_layout=True)
ax0 = ax[0]
ax1 = ax[1]
ax2 = ax[2]
ax3 = ax[3]

xmin = datetime.datetime(2010, 8, 25)
xmax = datetime.datetime(2022, 8, 7)

#mean speeds
ax0.plot(windMean_time, airSpeedMean, label='Dutch Harbor Airport')
ax0.plot(metMean['Wspd_kts'], label='Skan MetStation')
ax0.legend()
ax0.set_title('Mean Wind Speeds')
ax0.set_xlim(xmin, xmax)
ax0.set_ylabel('Speed (Knots)')

#max speeds
ax1.plot(windMax_time, airSpeedMax, label='Dutch Harbor Airport')
ax1.plot(metMax['Wspd_kts'], label='Skan MetStation')
ax1.legend()
ax1.set_title('Max Wind Speeds')
ax1.set_xlim(xmin, xmax)
ax1.set_ylabel('Speed (Knots)')

#mean direction
ax2.plot(windMean_time, airDirMean, label='Dutch Harbor Airport')
ax2.plot(metMean['Wdir_degT'], label='Skan MetStation')
ax2.legend()
ax2.set_title('Mean Wind Direction')
ax2.set_xlim(xmin, xmax)
ax2.set_ylabel('Direction (degT)')

#Direction at maximum wind speed
ax3.plot(windMax_time, airDirMax, label='Dutch Harbor Airport')
ax3.plot(metMax['Wdir_degT'], label='Skan MetStation')
ax3.legend()
ax3.set_title('Direction at Maximum Wind Speed')
ax3.set_xlim(xmin, xmax)
ax3.set_ylabel('Direction (degT)')

#%% Plot Dutch Harbor Wind Speed/Direction With Events Marked

figs,axs = plt.subplots(2,1, dpi=150, constrained_layout=True)
ax0=axs[0]
ax1=axs[1]
# ax2=axs[2]
# ax3=axs[3]

#### Dutch Harbor Wind 2014-2022
ax00 = ax0.twinx()
#direction
ax0.axhspan(300,360, color='k', alpha=0.1, label='300 degT',zorder=0)
p1 = ax0.plot(windMean_time, airDirMean, color='sandybrown',alpha=0.7,label='direction',zorder=1)
ax0.set_ylabel('Direction (degT)')
ax0.set(ylim=(0,360))

#speed
p2 = ax00.plot(windMax_time, airSpeedMax, color='tab:blue', linewidth=1.5, label='speed',zorder=1)
ax00.axvspan(datetime.date(2014,11,7), datetime.date(2014,11,9), color='k', alpha=0.2,zorder=3)
ax00.axvspan(datetime.date(2015,12,12), datetime.date(2015,12,14), color='k', alpha=0.2,zorder=3)
ax00.axvspan(datetime.date(2017,11,26), datetime.date(2017,11,28), color='k', alpha=0.2,zorder=3)
ax00.axvspan(datetime.date(2021,10,18), datetime.date(2021,10,24), color='k', alpha=0.2,zorder=3)
ax00.set_ylabel('Speed (knots)')
leg = p1 + p2
lab = [l.get_label() for l in leg]
ax0.legend(leg, lab, loc=0)
ax0.grid(color='lightgray', linestyle=':', linewidth=0.7, axis='x', which='major')
ax0.set_xlim(datetime.date(2014,1,1), datetime.date(2022,8,1))



#### Dutch Harbor vs Skan Wind (2021)
ax01 = ax1.twinx()
#speed
ax1.axhspan(300,360, color='k', alpha=0.1, label='300 degT',zorder=0)
p1 = ax1.plot(windMean_time, airDirMean, color='sandybrown', alpha=0.7,label='Dutch Direction',zorder=1)
p11 = ax1.plot(metMax['Wdir_degT'], color='sienna', label='Skan Direction', zorder=1)
ax1.set_ylabel('Direction (degT)')
ax1.set(ylim=(0,360))

#direction
p2 = ax01.plot(windMax_time, airSpeedMax, color='tab:blue', label='Dutch Speed',linewidth=2, zorder=2)
p22 = ax01.plot(metMax['Wspd_kts'], color='darkblue', label='Skan Speed', linewidth=2, zorder=2)
ax01.axvspan(datetime.date(2021,10,20), datetime.date(2021,10,23), color='k', alpha=0.2,zorder=3)
ax01.set_ylabel('Speed (knots)')
leg = p1 + p11 + p2 + p22
lab = [l.get_label() for l in leg]
ax1.legend(leg, lab, loc=0)
ax1.grid(color='lightgray', linestyle=':', linewidth=0.7, axis='x', which='major')
ax1.set_xlim(datetime.date(2021,8,1), datetime.date(2022,8,1))





#%% Compare Skan Buoy Data to NDBC Buoy

fig, ax = plt.subplots(4,1, figsize=(6,10), constrained_layout=True, dpi=150)
ax0 = ax[0]
ax1 = ax[1]
ax2 = ax[2]
ax3 = ax[3]

#mean wave heihgt
ax0.plot(ndbc_time, ndbcHeightMean, label= 'NDBC')
ax0.plot(waveMean_in['height'], label = 'Skan Inner Buoy')
ax0.plot(waveMean_out['height'], label='Skan Outer Buoy')
ax0.set_title('Mean Wave Heights')
ax0.legend()


#max wave height
ax1.plot(ndbc_time, ndbcHeightMax, label= 'NDBC')
ax1.plot(waveMax_in['height'], label = 'Skan Inner Buoy')
ax1.plot(waveMax_out['height'], label='Skan Outer Buoy')
ax1.set_title('Max Wave Heights')
ax1.legend()


#mean wave direction
ax2.plot(ndbc_time, ndbcDirMean, label= 'NDBC')
ax2.set_title('Mean Wave Dir')
ax2.legend()


#max wave direction
ax3.plot(ndbc_time, ndbcDirMax, label= 'NDBC')
ax3.set_title('Max Wave Dir')
ax3.legend()
