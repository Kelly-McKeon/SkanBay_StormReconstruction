#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 13:35:43 2024

@author: kellymckeon
"""

import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42
import pandas as pd
import time
import random

mac = r"/Users/kellymckeon/Library/CloudStorage/OneDrive-WoodsHoleOceanographicInstitution/Documents/WHOI/Alaska/Codes/Coherent_Codes"
pc = r"C:/Users/kelly/OneDrive - Woods Hole Oceanographic Institution/Documents/WHOI/Alaska/Codes/Coherent_Codes"

# CHANGE DIRECTORY TO WHICHEVER COMPUTER YOU ARE ON
os.chdir(mac)

start_code = time.time()

#%% Load Data

#### grainsize data 
JPCC = np.genfromtxt(r"JPC05/For_Publication/core_data/JPC05_Composite_GS.csv", delimiter=",", skip_header=1)
GS_depth = JPCC[:, 0] 
d10 = JPCC[:, 1]
d50 = JPCC[:, 2]
d90 = JPCC[:, 3]

#### age data
#choose age model slumps based on the threshold set in the next section
#if you have not run the age model yet, comment out this section and return to it after running the section below
ages = np.genfromtxt(r"JPC05/For_Publication/core_data/JPC05_Composite_ages.txt", skip_header=1) 
ages_depth = ages[:,0]
ind = np.isin(ages_depth, GS_depth) #get rid of duplicate ages from slumps

# ages are in YBP
ages_med = ages[ind,3] 
ages_mean = ages[ind,4] 
ages_min = ages[ind,1] 
ages_max = ages[ind,2] 




#%% Visualize Grainsize Data

#### plot d50 histogram to get lay of the land
plt.figure()
plt.hist(d50, bins=20, histtype='bar', ec='k')
plt.xlabel('grainsize ($\mu$m)')
plt.ylabel('counts')
plt.title('JPCC d50')
plt.show()

#save the histogram output 
count, bins_count = np.histogram(d50, bins=20)


#### plot d50 PDF and CDF

#calculate PDF
pdf = count / sum(count)
#calculate CDF
cdf = np.cumsum(pdf)

#plot d50 distribution
plt.figure(dpi=150)
plt.plot(bins_count[1:], pdf, color="red", label="PDF")
plt.plot(bins_count[1:], cdf, label="CDF")
plt.axvline(np.std(d50) + np.median(d50), color='k', linestyle='dashed', label='d50 of 1$\sigma$')
plt.axvline(2*np.std(d50) + np.median(d50), color='gray', linestyle='dashed', label='d50 of 2$\sigma$')
plt.axvline(63, color='k', label='63 $\mu$m')
plt.ylabel('Distribution Function')
plt.xlabel('d50 ($\mu$m)')
plt.title('d50 Distribution')
plt.legend()
plt.show()


#%% Create Slump Threshold
sd50 = np.std(d50) #d50 standard deviation
med50 = np.median(d50) #d50 median
avg50 =np.mean(d50) #d50 mean
thr50 = 63 #use threshold for sand

# get points of indices greater than the threshold
tind = np.argwhere(d50 >= thr50) #save indices greater than the threshold
d50_2 = np.copy(d50) #make a new array because python back changes your array if you save it as the same name
d50_2[tind] = med50 #replace values above the threshold with the timeseries median


#### Get Event Depths to slump in Age Model
event_depths = GS_depth[tind]
event_gs = d50[tind]


#### Plot d50 with Slumps Marked (Sanity Check)
fig,ax = plt.subplots(figsize=(4,8),dpi=150)
ax.plot(d50, GS_depth, label = 'd50') #grainsize
ax.scatter(event_gs, event_depths, marker = '*', color='k', label='slumps') #slumps
ax.plot(d50_2, GS_depth, color='k', label= "d50 < 63$\mu$m") #grainsize below threshold
plt.axvline(x=thr50, label=str(int(thr50)),color='k', linestyle = 'dashed') #threshold
ax.invert_yaxis()
ax.set_ylabel('depth (cm)')
ax.set_xlabel('d50 ($\mu$m)')
ax.set_title('Slump Depths for Age Model')
ax.legend()
plt.show()

#### Save Slump Starts and Slump Ends
# open these variables and type them into the age model program
iterator = np.arange(np.size(event_depths))
starts = np.zeros(np.size(iterator)) * np.nan
ends = np.zeros(np.size(iterator)) * np.nan
slump_start = []
slump_end = []

#get cm difference between event depths
for i in iterator[:-1]:
    starts[i] = event_depths[i+1] - event_depths[i] 

for i in iterator:
    if i == 0: #skip first deposit
        pass
    else: 
        #if distance between this event and the last event is greater than 1, it's a new slump start
        if starts[i] == 1 and starts[i-1] > 1:
           slump_start.append(event_depths[i])
        #if distance between this event and the next event is greater than 1, it's a new slump end
        if starts[i] == 1 and starts[i+1] > 1:
            slump_end.append(event_depths[i+1])
        #the core ends with an event
        #in this case just append the last event depth to the end 
        if i == len(event_depths)-1:
            slump_end.append(event_depths[i])
            
###############################################################################           

##### Now go in rPlum and add slumps at all the event depths that are over the threshold for more than 1 cm

###############################################################################

#%% Make Grainsize Anomaly Timeseries

# 25cm window is centennial avg assuming 4 yr/cm accumulation rate
window_size = 25 

# Convert array of integers to pandas series because rolling windows is easier in Pandas
d50_series = pd.Series(d50_2) #use the timeseries with the outliers removed so they do not skew the moving median
  
# Get the windows
windows = d50_series.rolling(window_size)
  
# Create a series of moving averages of each window
mov = windows.mean()

#fill the nans left out of the moving average with the ts mean
moving_averages = mov.fillna(avg50)
  
# Convert pandas series back to numpy 
d50_mov = moving_averages.to_numpy()

#### Get Anomaly
anomaly50 = d50 - d50_mov #use the full timeseries not the one with outliers removed

#### Grainsize Anomaly (Sanity Check)
fig, ax = plt.subplots(1,2, figsize=(6,8),dpi=150,constrained_layout=True)
ax0=ax[0]
ax1=ax[1]

# d50 with anomaly
ax0.plot(d50, GS_depth, label = 'JPCC d50') #grainsize
ax0.plot(d50_mov, GS_depth, label = '25pt mov med') #100yr moving mean
ax0.invert_yaxis()
ax0.set_ylabel('depth (cm)')
ax0.set_xlabel('d50 ($\mu$m)')
ax0.set_title('JPCC d50')
ax0.legend()

# anomaly with potential event thresholds marked
ax1.plot(anomaly50, GS_depth) #anomaly
ax1.axvline(0, color='k', label='zero') #plot zero line
ax1.axvline(np.median(anomaly50)+np.std(anomaly50), color='k', linestyle = 'dashed', label='1$\sigma$')
ax1.axvline(np.median(anomaly50)+2*np.std(anomaly50), color='k', linestyle = 'dotted', label='2$\sigma$')
ax1.invert_yaxis()
ax1.set_ylabel('depth (cm)')
ax1.set_xlabel('d50 Anomaly ($\mu$m)')
ax1.set_title('JPCC d50 Anomaly')
ax1.legend()


#%% Create Event Threshold

#### Visualize Anomalies

#increase bin size until the mode is centered on zero
plt.figure(dpi=150)
plt.hist(anomaly50,bins=50,histtype='bar',ec='k')
plt.title('d50 anomaly ($\mu$m)')
plt.show()

#save the output
count, bins_count = np.histogram(anomaly50, bins=50)

#make PDF
pdf = count / (sum(count) * np.diff(bins_count))
#make CDF
cdf = np.cumsum(pdf)

#### Get Anomaly Summary Statistics
sd50a = np.std(anomaly50)
med50a = np.median(anomaly50)
avg50a =np.mean(anomaly50)
uppera = np.percentile(anomaly50, 90)
lowera = np.percentile(anomaly50, 83) #vary this one to choose event threshold
athresh = anomaly50 >= lowera

#choose 85 as the event threshold because it is at the bend point of the PDF plot

#### Plot PDF and CDF
plt.figure(dpi=150)
plt.plot(bins_count[1:], pdf, color="red", label="PDF")
plt.axvline(lowera, color='k', linestyle = 'dashed', label = 'gs of 83rd')
plt.axvline(np.percentile(anomaly50, 85), color = 'k', label = 'gs of 85th')
plt.axvline(uppera, color='gray', linestyle = 'dashed', label = 'gs of 90th')
plt.plot(bins_count[1:], cdf, label="CDF")
plt.xlabel('d50 anomaly ($\mu$m)')
plt.ylabel('distribution functions')
plt.title('25pt filter anomaly')
plt.legend()
plt.show()


#### Plot Anomaly with Threshold Options Marked 
fig, ax = plt.subplots(1,2,figsize=(4,8),constrained_layout=True,dpi=150)
ax0=ax[0]
ax1=ax[1]

ax0.plot(d50, GS_depth, label='d50')#grainsize plotted against depth
ax0.plot(d50_mov, GS_depth, label = 'd50 filter') #moving median plotted with depth
ax0.scatter(event_gs, event_depths, marker = '*', color='k')
ax0.invert_yaxis()
ax0.set_ylabel('depth (cm)')
ax0.set_xlabel('d50 ($\mu$m)')

ax1.plot(anomaly50, ages_med, label = 'anomaly') #anomaly plotted against age
ax1.scatter(anomaly50[athresh], ages_med[athresh], marker = '*', color='k') #depths above the event threshold
ax1.invert_yaxis()
ax1.axvline(lowera, color = 'k', linestyle = 'dashed', linewidth=1, label= '85th percentile anomaly')
ax1.set_xlabel('coarse anomaly ($\mu$m)')
ax1.set_ylabel('median age (YBP)')
plt.suptitle('d50 Events')
plt.xlim(0,300)
plt.show()


#### Plot Histogram with PDF over it
plt.figure(dpi=150)
plt.plot(bins_count[1:], pdf, color='red', label='PDF')
plt.axvline(lowera, color='k', label='85th percentile')
plt.hist(anomaly50,bins=70, density=True, histtype='bar',ec='k')
plt.rc('axes', axisbelow=True)
plt.grid(color='lightgray', linestyle=':', linewidth=0.7, axis='both', which='major')
plt.xlabel('d50 anomaly ($\mu$m)')
plt.xlim(-30,100)
plt.ylabel('density')
plt.legend()

#%% Make Boolean Event Index

new_array = np.ones(np.size(anomaly50)) #make array that will populate with Boolean index

event_depth_array = [] #make empty depth array to save event depths in
event_gs_array = [] #make empty gs array to save event grainsizes in
event_anom_array = [] #make empty anomaly array
event_time_array = [] #make empty time array

iterator = np.arange(0,np.size(anomaly50))

#### Set Each Depth as Event or Not
#use the lower threshold defined above

for i in iterator:
    if anomaly50[i] < lowera:
        new_array[i] = 0
        
    #only code it as an event if the cm before was NOT an event
    elif anomaly50[i] >= lowera and anomaly50[i-1] < lowera: 
        new_array[i] = 1
        event_depth_array.append(GS_depth[i])
        event_gs_array.append(d50[i])
        event_anom_array.append(anomaly50[i])
        event_time_array.append(ages_med[i])
   
    #if both consecutive cm are events, only code the first one as an event
    elif anomaly50[i] >= lowera and anomaly50[i-1] >= lowera: 
        new_array[i] = 0
        
    bool_anom = new_array #make boolean event array
        
#%% Remove Tsunamis from Event Depths 

del event_depth_array[19] #deposit at 518cm
del event_depth_array[55] #deposit at 1415cm

del event_gs_array[19] #deposit at 518cm
del event_gs_array[55] #deposit at 1415cm

del event_anom_array[19] #deposit at 518cm
del event_anom_array[55] #deposit at 1415cm

del event_time_array[19] #deposit at 518cm
del event_time_array[55] #deposit at 1415cm


#set tsunamis to non-events
bool_anom[517] = 0 
bool_anom[1415] = 0

#%% Plot Events

#### Plot d50 with discrete events marked
fig, ax = plt.subplots(1,2,figsize=(8,8),constrained_layout=True,dpi=150)
ax0=ax[0]
ax1=ax[1]

#plot raw grainsize with events marked as sanity check (should be just the first cm of event)
ax0.plot(d50, GS_depth, label='d50')
ax0.invert_yaxis()
ax0.set_ylabel('depth (cm)')
ax0.set_xlabel('d50 ($\mu$m)')
ax0.scatter(event_gs_array, event_depth_array, marker='*', color='k', label='events')
ax0.axvline(np.min(event_gs_array), color = 'k', linestyle = 'dashed', label='minimum event GS')
ax0.legend(bbox_to_anchor=(0.6, 0.95))

#plot grainsize anomaly with events marked at median ages
ax1.plot(anomaly50, ages_med, label = 'anomaly')
ax1.scatter(event_anom_array, event_time_array, marker='*', color='k', label='events')
ax1.invert_yaxis()
ax1.axvline(lowera, color = 'k', linestyle = 'dashed', label= '85th prctile')
ax1.set_xlabel('coarse anomaly ($\mu$m)')
ax1.set_ylabel('median age (YBP)')
ax1.legend(bbox_to_anchor=(0.6, 0.95))
plt.xlim(0,300)
plt.show()

#### plot grainsize anomaly with event cutoff
fig, ax = plt.subplots(1,1, dpi=150)
ax.plot(ages_med, anomaly50, label = 'anomaly', zorder=1)
ax.scatter(event_time_array, event_anom_array, marker='*', color='k', label='events', zorder=2)
ax.axhline(lowera, color = 'k', linestyle = 'dashed', linewidth=1, label= '85th prctile',zorder=3)
ax.set_ylabel('coarse anomaly ($\mu$m)')
ax.set_xlabel('median age (YBP)')
ax.legend()
plt.show()

#%% Plot Event Threshold with Raw d50

mov_anomaly = d50_mov + lowera

fig,ax = plt.subplots(1,1, figsize=(8,4), dpi=150)
ax.plot(ages_med, d50, label='d50', zorder=1)
ax.plot(ages_med[25:], mov_anomaly[25:], linestyle='--', linewidth=1, color='k',label='threshold', zorder=2)
ax.scatter(event_time_array, event_gs_array, marker='*', color='k', label='events', zorder=3)
ax.set_xlabel('median age (YBP)')
ax.set_ylabel('d50 ($\mu$m)')
plt.grid(color='lightgray', linestyle=':', linewidth=0.7, axis='both', which='major')
ax.legend()


#%% Save Event Depths For Export to R
np.savetxt(r"JPC05/For_Publication/core_data/event_depth.csv", event_depth_array, delimiter=',')
np.savetxt(r"JPC05/For_Publication/core_data/event_gs.csv", event_gs_array, delimiter=',')

#once you have event depths go to rPlum and export age distributions for each of these depths

#%% Load Age Distributions from Bacon
#ages are in YBP
age_dist = np.genfromtxt(r"JPC05/For_Publication/core_data/JPC05_Age_Distributions.csv", delimiter=",", skip_header=1)
age_dist = age_dist[:,1:] #remove first index column of nans


#%% Make Dataset Annual

eventyrs=[]
#populate new list with event years based on median age
for i in range(0, len(bool_anom)):
    if bool_anom[i] == True:
        eventyrs.append(ages_med[i]) 
        
eventyrs = np.array(eventyrs) #convert to array
eventyrs = eventyrs.astype(int) #round to nearest year (no decimals)

#make new array with evently spaced years
allyrs = np.arange(min(ages_med), max(ages_med), 1, dtype=int) #boundaries are the minimum and maximum median ages, stepping by 1 year at a time
strmyrs = np.isin(allyrs, eventyrs) #this creates a binary index of the same size as allyrs, where if the year in allyrs is a storm year (in yrs) it's coded as 1

#### Stem Plot of Storm Years
plt.figure(figsize=(6,4), dpi=150)
plt.stem(allyrs,strmyrs, 'k', markerfmt =" ", basefmt='k-')
plt.xlabel('years before present')
plt.ylabel('storm')
plt.xlim(-100,3000) #cut it off at 3000ybp because of the sea level transgression

#save event years
np.savetxt(r"JPC05/For_Publication/core_data/eventyrs.csv", eventyrs, delimiter=',')
# save storm years (binary)
np.savetxt(r"JPC05/For_Publication/core_data/storm_years.csv", strmyrs, delimiter=',')
# save allyrs
np.savetxt(r"JPC05/For_Publication/core_data/all_years.csv", allyrs, delimiter=',')

#%% Gaussian Kernel Rate Estimation

start_time = time.time()
#set bandwidth
h = 50 #500 year bandwidth = 100 year window

#### Create Function
#create a function for easy bootstrapping later

def rate_parameter(allyrs, eventyrs, h):
    #initialize the rate parameter
    l = np.nan*(np.zeros(len(allyrs))) #l is lambda as in Mudelsee Eq 3.16
    
    for i in range(0,len(allyrs)): #loop through every year of the dataset
        l[i] = 0 
        for j in range(0, len(eventyrs)):
            y = (allyrs[i] - eventyrs[j]) / h
            l[i] = l[i] + (np.power((2*np.pi), -0.5) * np.exp((-y**2)/2))/h
    return(l*100) # the rate parameter is the # of events per year, so to get storms/century multiply by 100

#### Calculate Rate Parameters
l25 = rate_parameter(allyrs, eventyrs, 25)
l50 = rate_parameter(allyrs, eventyrs, 50)
l100 = rate_parameter(allyrs, eventyrs, 100)

print('time to run kernel rate estimation is ', (time.time() - start_time)/60, ' minutes')  


#%% Plot Median Age Rate Parameter with MCA and LIA Marked

plt.figure(dpi=150)
plt.stem(allyrs,strmyrs, 'k', markerfmt =" ", basefmt=" ")
plt.plot(allyrs, l50,'tab:blue', linewidth=2)
plt.xlabel('Years Before Present (YBP)')
plt.ylabel('Storm Frequency ($\lambda$)')
plt.axvspan(xmin=700, xmax=1000, color='tab:red', alpha=0.1, label='MCA')
plt.axvspan(xmin=100, xmax=650, color='tab:blue', alpha=0.1, label='LIA')
plt.xlim(-100,3000)
plt.legend()

#%% Uniform Kernel Rate Estimation

window=50 #window of 50 is 100 year frequency


#### Define Uniform Function
def uniform_kernel(allyrs, strmyrs, window):
    ages_temp = allyrs
    new = np.nan*np.ones(len(allyrs)) #new array to populate              
    for i in range(0,len(ages_temp)):
        if ages_temp[i] >= min(ages_temp) + window and ages_temp[i] <= max(ages_temp) - window: #age needs to not be within 50 years of the timeseries edge
            i2 = np.argwhere(np.logical_and(ages_temp >= (ages_temp[i] - window), ages_temp <= (ages_temp[i]+window))) #get all the indices +/- 50 years from each cm
            freq = sum(strmyrs[i2]) #test if there was a storm in each cm associated with the 100 year age range
            new[i] = freq
        elif ages_temp[i] < min(ages_temp) + window:
            new[i] = 0
        elif ages_temp[i] > max(ages_temp)-window:
            new[i] = 0
            
    return(new)
            
#### Get Uniform Frequency
freq_array_100 = uniform_kernel(allyrs, strmyrs, window) #make freq array 

#%% Save Median Ages Storm Frequencies
ages_CE = 1950-(allyrs)
df = pd.DataFrame(ages_CE)
df.to_csv(r"JPC05/For_Publication/core_data/age_CE.csv", index=False) #CE ages

df = pd.DataFrame(freq_array_100)
df.to_csv(r"JPC05/For_Publication/core_data/unif_freq.csv", index=False) #uniform frequency

df = pd.DataFrame(l50)
df.to_csv(r"JPC05/For_Publication/core_data/med_age_freq.csv", index=False) #Gaussian frequency

#%% Plot Comparsion of Frequency Calculation Methods 

plt.figure(figsize=(8,4), dpi=300)
plt.plot(allyrs, freq_array_100, label='Uniform Kernel', color='k', linewidth=2)
plt.plot(allyrs, l25, label='h = 25', color= 'tab:blue')
plt.plot(allyrs, l50, label='h = 50', color= 'tab:orange')
plt.plot(allyrs, l100, label='h = 100', color= 'tab:green')
plt.xlabel('YBP')
plt.ylabel('rate parameter')
plt.legend()
plt.title('kernel comparison')


#%% Load in Bootstrap Data if not Bootstrapping
l_matrix_noages = np.genfromtxt(r"JPC05/For_Publication/core_data/l_matrix.csv", delimiter=",", skip_header=0)
l_matrix_ages = np.genfromtxt(r"JPC05/For_Publication/core_data/l_matrix_ages.csv", delimiter=",", skip_header=0) 
lp_matrix_noages = np.genfromtxt(r"JPC05/For_Publication/core_data/lp_matrix.csv", delimiter=",", skip_header=0)
lp_matrix_ages = np.genfromtxt(r"JPC05/For_Publication/core_data/lp_matrix_ages.csv", delimiter=",", skip_header=0) 

#%% Bootstrap Gaussian Kernel Rate Estimation without Age Uncertainty (~30 min)

start_time = time.time()

#set number of iterations
iterations = 1000
#set window size
h = 50
#initialize empty matrix
l_matrix_noages = np.zeros((iterations, len(allyrs)))

for ii in range(0,iterations):
  
    #create a new set of event years with replacement from event years dataset
    rmax = len(eventyrs) #maximum index
    rmin = 0 #minimum index
    event_ind = np.random.randint(rmin, high=rmax, size=len(eventyrs)) #rmax is exclusive of the final element 
    eventyrs_boot = eventyrs[event_ind]
    
    #perform rate parameter analysis on the new timeseries of event years
    l_temp = rate_parameter(allyrs, eventyrs_boot, h)
    
    l_matrix_noages[ii,:] = l_temp
    
print('time to run bootstrapping is ', (time.time() - start_time)/60, ' minutes')   

#%% Plot Rate Parameter without Age Uncertainties  

#### Get Error Bars
mean_rate_noages = np.median(l_matrix_noages, axis=0)
low_rate2_noages = np.percentile(l_matrix_noages, 5, axis=0)
low_rate1_noages = np.percentile(l_matrix_noages, 16, axis=0)
high_rate2_noages = np.percentile(l_matrix_noages, 95, axis=0)
high_rate1_noages = np.percentile(l_matrix_noages, 84, axis=0)

#### Plot 
fig,ax = plt.subplots(figsize=(12,4), dpi=150)
# ax.plot(allyrs, freq_array_100, label='Uniform Kernel', color='navy')
ax.plot(allyrs, mean_rate_noages, color='k', label = 'bootstrapped mean')
ax.fill_between(allyrs, low_rate2_noages, high_rate2_noages, color = 'k', alpha=0.3, label='+/- 2sd')
ax.fill_between(allyrs, low_rate1_noages, high_rate1_noages, color='k', alpha=0.5, label = '+/- 1sd')
ax.plot(allyrs, l50, label='median event age')
ax.legend()
ax.set_xlabel('Years Before Present (YBP)')
ax.set_ylabel('Storm Frequency ($\lambda$)')
ax.set_title('Boostrapping without Age Uncertainty')


#### Get Lambda Statistics

#rate parameter max
lmax_noages = np.max(l_matrix_noages, axis=1)
lmax_mean_noages = np.mean(lmax_noages)
lmax_5_noages = np.percentile(lmax_noages, 5)
lmax_95_noages = np.percentile(lmax_noages, 95)
lmax_std_noages = np.std(lmax_noages)


#year of maximum rate parameter
lmax_ind_array_noages = np.argmax(l_matrix_noages, axis=1)
year_max_noages=np.zeros(len(lmax_ind_array_noages))
for i in range(0,len(lmax_ind_array_noages)):
    year_max_noages[i] = allyrs[lmax_ind_array_noages[i]]

yrmax_mean_noages = np.mean(year_max_noages)
yrmax_5_noages = np.percentile(year_max_noages, 5)
yrmax_95_noages = np.percentile(year_max_noages, 95)
yrmax_std_noages = np.std(year_max_noages)

 

#%% Bootstrap Gaussian Kernel Rate Estimation with Age Uncertainties (~30 min)


start_time = time.time()

#set number of iterations
iterations = 1000
#set window size
h = 50
#initialize empty matrix
l_matrix_ages = np.zeros((iterations, len(allyrs)))
event_matrix_bool = np.nan * np.zeros((iterations,len(allyrs)))

for ii in range(0,iterations):
    
    #draw a random event age from the age probability distributions
    amax = len(age_dist) #maximum index
    amin = 0 #minimum index
    age_ind = np.random.randint(amin, high=amax, size=1) #get random index
    
    #save alternative set of event years from the random age model iteration
    eventyrs_temp = age_dist[age_ind,:]
    eventyrs_int = np.squeeze(eventyrs_temp.astype(int)) #convert to integer
    
    #create a new set of event years with replacement from the intermediate event years 
    #keeping the same total number of events
    #preserving stratigraphic order
    rmax = len(eventyrs) #maximum index
    rmin = 0 #minimum index
    event_ind = np.random.randint(rmin, high=rmax, size=len(eventyrs)) #rmax is exclusive of the final element so add 1
    eventyrs_boot = eventyrs_int[event_ind]
    
    strmyrs_boot = np.isin(allyrs, eventyrs_boot) #boolean index of event years
    
    #perform rate parameter analysis on the new timeseries of event years
    l_temp = rate_parameter(allyrs, eventyrs_boot, h)
    
    l_matrix_ages[ii,:] = l_temp
    event_matrix_bool[ii,:] = strmyrs_boot

print('time to run bootstrapping is ', (time.time() - start_time)/60, ' minutes')    

#%% Get Lambda Statistics

#rate parameter max
lmax = np.max(l_matrix_ages, axis=1)
lmax_mean = np.mean(lmax)
lmax_5 = np.percentile(lmax, 5)
lmax_95 = np.percentile(lmax, 95)


#year of maximum rate parameter
lmax_ind_array = np.argmax(l_matrix_ages, axis=1)
year_max=np.zeros(len(lmax_ind_array))
for i in range(0,len(lmax_ind_array)):
    year_max[i] = allyrs[lmax_ind_array[i]]

yrmax_mean = np.mean(year_max)
yrmax_5 = np.percentile(year_max, 5)
yrmax_95 = np.percentile(year_max, 95)
yrmax_std = np.std(year_max)


#%% Plot Rate Parameter with Age Uncertainty 

#### Get Error Bars
mean_rate_ages = np.median(l_matrix_ages, axis=0)
low_rate2_ages = np.percentile(l_matrix_ages, 5, axis=0)
low_rate1_ages = np.percentile(l_matrix_ages, 16, axis=0)
high_rate2_ages = np.percentile(l_matrix_ages, 95, axis=0)
high_rate1_ages = np.percentile(l_matrix_ages, 84, axis=0)

#### Plot 
fig,ax = plt.subplots(figsize=(12,4), dpi=150)
ax.plot(allyrs, mean_rate_ages, color='k', label = 'bootstrapped mean')
ax.fill_between(allyrs, low_rate2_ages, high_rate2_ages, color = 'k', alpha=0.3, label='2$\sigma$')
ax.fill_between(allyrs, low_rate1_ages, high_rate1_ages, color='k', alpha=0.5, label = '1$\sigma$')
ax.plot(allyrs, freq_array_100, label='Uniform Kernel', color='navy')
ax.plot(allyrs, l50, label='median event age')
ax.legend()
ax.set_xlabel('Years Before Present (YBP)')
ax.set_ylabel('Storm Frequency ($\lambda$)')
ax.set_title('Boostrapping with Age Uncertainty')
ax.set_xlim(-150,8150)


#### Get Stats

rmean = np.mean(mean_rate_ages[0:3000]) #record mean
rstd = np.std(mean_rate_ages[0:3000]) #record std

premax_mean = np.mean(mean_rate_ages[0:870+72]) #premax mean
postmax_mean = np.mean(mean_rate_ages[870+72:3000]) #postmax mean
first100_mean = np.mean(mean_rate_ages[0:100]) #last century mean

#%% Make Active Intervals

#### Make Thresholds

#without ages
high_thresh_meanrate_noages = np.percentile(mean_rate_noages[0:3000], 84) #1sigma
low_thresh_meanrate_noages = np.percentile(mean_rate_noages[0:2800], 16) #cut low thresh where rate goes to zero

#ages
high_thresh_meanrate_ages = np.percentile(mean_rate_ages[0:3000], 84) #1sigma
low_thresh_meanrate_ages = np.percentile(mean_rate_ages[0:2800], 16)#cut low thresh where rate goes to zero

#### Analysis C
ind = np.argwhere(high_rate2_noages >= high_thresh_meanrate_noages) #where 90th percentile lambda is above 1sigma of all lambdas
activeyrs = allyrs[ind]
active_starts = []
active_ends = []

active_starts.append(activeyrs[0])

for i in range(0,len(activeyrs)-1): #already added the last point to the ends so don't go to the end here or there will be an index problem
    if activeyrs[i+1] - activeyrs[i] == 1:
        pass
    elif activeyrs[i+1] - activeyrs[i] > 1:
        active_starts.append(activeyrs[i+1]) #append only works with lists 
        active_ends.append(activeyrs[i])

active_ends.append(activeyrs[-1]) 

#save active years without ages
activeyrs_noages = activeyrs #active years
active_starts_noages = np.array(active_starts) #beginning of active intervals
active_ends_noages = np.array(active_ends) #end of active intervals


#### Analysis D
ind = np.argwhere(high_rate2_ages >= high_thresh_meanrate_ages)  #where 90th percentile lambda is above 1sigma of all lambdas
activeyrs = allyrs[ind]
active_starts = []
active_ends = []

active_starts.append(activeyrs[0])

for i in range(0,len(activeyrs)-1): #already added the last point to the ends so don't go to the end here or there will be an index problem
    if activeyrs[i+1] - activeyrs[i] == 1:
        pass
    elif activeyrs[i+1] - activeyrs[i] > 1:
        active_starts.append(activeyrs[i+1]) #append only works with lists 
        active_ends.append(activeyrs[i])

active_ends.append(activeyrs[-1]) 

activeyrs_ages = activeyrs
active_starts_ages = np.array(active_starts)
active_ends_ages = np.array(active_ends)

#### Analysis B
ind = np.argwhere(mean_rate_noages >= high_thresh_meanrate_noages)
activeyrs = allyrs[ind]
active_starts = []
active_ends = []

active_starts.append(activeyrs[0])

for i in range(0,len(activeyrs)-1): #already added the last point to the ends so don't go to the end here or there will be an index problem
    if activeyrs[i+1] - activeyrs[i] == 1:
        pass
    elif activeyrs[i+1] - activeyrs[i] > 1:
        active_starts.append(activeyrs[i+1]) #append only works with lists 
        active_ends.append(activeyrs[i])

active_ends.append(activeyrs[-1]) 

activeyrs_mean = activeyrs
active_starts_mean = np.array(active_starts)
active_ends_mean = np.array(active_ends)



#%% Make Quiet Intervals

#### Analysis C
ind = np.argwhere(low_rate2_noages[0:2800] <= low_thresh_meanrate_noages)  #where 90th percentile lambda is below 1sigma of all lambdas
quietyrs = allyrs[ind]
quiet_starts = []
quiet_ends = []

quiet_starts.append(quietyrs[0])

for i in range(0,len(quietyrs)-1): #already added the last point to the ends so don't go to the end here or there will be an index problem
    if quietyrs[i+1] - quietyrs[i] == 1:
        pass
    elif quietyrs[i+1] - quietyrs[i] > 1:
        quiet_starts.append(quietyrs[i+1]) #append only works with lists 
        quiet_ends.append(quietyrs[i])

quiet_ends.append(quietyrs[-1])   

quietyrs_noages = quietyrs
quiet_starts_noages = np.array(quiet_starts)
quiet_ends_noages = np.array(quiet_ends)


#### Analysis B
ind = np.argwhere(low_rate2_ages[0:2800] <= low_thresh_meanrate_ages)
quietyrs = allyrs[ind]
quiet_starts = []
quiet_ends = []

quiet_starts.append(quietyrs[0])

for i in range(0,len(quietyrs)-1): #already added the last point to the ends so don't go to the end here or there will be an index problem
    if quietyrs[i+1] - quietyrs[i] == 1:
        pass
    elif quietyrs[i+1] - quietyrs[i] > 1:
        quiet_starts.append(quietyrs[i+1]) #append only works with lists 
        quiet_ends.append(quietyrs[i])

quiet_ends.append(quietyrs[-1])   

quietyrs_ages = quietyrs
quiet_starts_ages = np.array(quiet_starts)
quiet_ends_ages = np.array(quiet_ends)


#### Analysis C
#mean rate no ages is equivalent to median age timeseries
ind = np.argwhere(mean_rate_noages[0:2800] <= low_thresh_meanrate_noages)
quietyrs = allyrs[ind]
quiet_starts = []
quiet_ends = []

quiet_starts.append(quietyrs[0])

for i in range(0,len(quietyrs)-1): #already added the last point to the ends so don't go to the end here or there will be an index problem
    if quietyrs[i+1] - quietyrs[i] == 1:
        pass
    elif quietyrs[i+1] - quietyrs[i] > 1:
        quiet_starts.append(quietyrs[i+1]) #append only works with lists 
        quiet_ends.append(quietyrs[i])

quiet_ends.append(quietyrs[-1])   

quietyrs_mean = quietyrs
quiet_starts_mean = np.array(quiet_starts)
quiet_ends_mean = np.array(quiet_ends)

#%%  Plot Active/Quiet Intervals

#plot only the rate maximum and minimum 
#hard code the one active interval to display so the plot doesn't get too messy
#get the values from the starts and ends variables made above 

plt.figure(figsize=(8,4), dpi=150)
plt.axvspan(xmin=758, xmax=929, color='darkolivegreen', alpha=0.3, label='Analysis B Active Interval',zorder=2)
plt.axvspan(xmin=515, xmax=1055, color='darkolivegreen', alpha=0.2, label='Analysic C 1.65$\sigma$ Active Interval',zorder=2)
plt.axvspan(xmin=188, xmax=1103, color='darkolivegreen', alpha=0.1, label='Analysis D 1.65$\sigma$ Active Interval',zorder=2)
plt.axvspan(xmin=1876, xmax=2038, color='firebrick', alpha=0.3, label='Analysis B Quiet Interval',zorder=2)
plt.axvspan(xmin=1822, xmax=2727, color='firebrick', alpha=0.2, label='Analysis C 1.65$\sigma$ Quiet Interval',zorder=2)
plt.axvspan(xmin=1452, xmax=2727, color='firebrick', alpha=0.1, label='Analysis D 1.65$\sigma$ Quiet Interval',zorder=2)
plt.axhline(high_thresh_meanrate_noages, color='dimgray', linestyle='dashed', label='Analysis C Active Threshold', zorder=3)
plt.axhline(high_thresh_meanrate_ages, color='dimgray', label='Analysis D Active Threshold', zorder=3)
plt.plot(allyrs, mean_rate_ages, color='k', label = 'Analysis D Mean',zorder=3)
plt.plot(allyrs, l50, label='Analysis B Mean', color='k', linestyle='dashed',zorder=3)
plt.fill_between(allyrs, low_rate2_ages, high_rate2_ages, color = 'lightgray', edgecolor='darkgray', linewidth=0.5, label='Analysis D 1.65$\sigma$',zorder=1)
plt.fill_between(allyrs, low_rate1_ages, high_rate1_ages, color = 'darkgray', edgecolor='gray',  label = 'Analysis D 1$\sigma$',zorder=1)
plt.grid(color='lightgray', linestyle=':', linewidth=0.7, axis='both', which='major', zorder=0)
plt.xlabel('Years Before Present (YBP)')
plt.ylabel('Storm Frequency ($\lambda$)')
plt.xlim(-100, 3000)
# plt.gca().invert_xaxis() 
plt.legend()
plt.title('Active/Inactive Intervals for Bootstrapped Ages')


#%% Plot Analysis B,C,D (3000 years)
fig, ax = plt.subplots(2,1, dpi=150)
ax0 = ax[0]
ax1 = ax[1]

ax0.stem(allyrs, strmyrs, 'k', markerfmt =" ", basefmt=" ")
ax0.plot(allyrs, l50, linestyle='--', color='k', linewidth=2, label='Analysis B', zorder=2)
ax0.plot(allyrs, mean_rate_ages, color='k', label = 'Analysis D Mean', zorder=2)
ax0.fill_between(allyrs, low_rate2_ages, high_rate2_ages, color = 'lightgray', edgecolor='darkgray', linewidth=0.5, label='Analysis D 1.65$\sigma$', zorder=1)
ax0.fill_between(allyrs, low_rate1_ages, high_rate1_ages, color='darkgray', edgecolor='gray', linewidth=0.5, label = 'Analysis D 1$\sigma$', zorder=1)
ax0.set_xlabel('Years Before Present (YBP)')
ax0.set_ylabel('Storm Frequency ($\lambda$)')
ax0.axvspan(xmin=700, xmax=1000, color='tab:red', alpha=0.1, label='MCA', zorder=0)
ax0.axvspan(xmin=100, xmax=650, color='tab:blue', alpha=0.1, label='LIA', zorder=0)
ax0.grid(color='lightgray', linestyle=':', linewidth=0.7, axis='both', which='major')
ax0.set_xlim(-100,3000)
ax0.legend()


ax1.stem(allyrs, strmyrs, 'k', markerfmt =" ", basefmt=" ")
ax1.plot(allyrs, l50, linestyle='--', color='k', linewidth=2, label='Analysis B', zorder=2)
ax1.plot(allyrs, mean_rate_noages, color='k', label = 'Analysis C Mean', zorder=2)
ax1.fill_between(allyrs, low_rate2_noages, high_rate2_noages, color = 'lightgray', edgecolor='darkgray', linewidth=0.5, label='Analysis C 1.65$\sigma$', zorder=1)
ax1.fill_between(allyrs, low_rate1_noages, high_rate1_noages, color='darkgray', edgecolor='gray', linewidth=0.5, label = 'Analysis C 1$\sigma$', zorder=1)
ax1.set_xlabel('Years Before Present (YBP)')
ax1.set_ylabel('Storm Frequency ($\lambda$)')
ax1.axvspan(xmin=700, xmax=1000, color='tab:red', alpha=0.1, label='MCA', zorder=0)
ax1.axvspan(xmin=100, xmax=650, color='tab:blue', alpha=0.1, label='LIA', zorder=0)
ax1.grid(color='lightgray', linestyle=':', linewidth=0.7, axis='both', which='major')
ax1.set_xlim(-100,3000)
ax1.legend()


#%% Plot Analysis A-D (8000 Years)
fig, ax = plt.subplots(2,1, dpi=150)
ax0 = ax[0]
ax1 = ax[1]

# ax0.stem(allyrs, strmyrs, 'k', markerfmt =" ", basefmt=" ")
ax0.plot(allyrs, freq_array_100, label='Analysis A', color='navy',zorder=2)
ax0.plot(allyrs, l50, color='tab:blue', linewidth=2, label='Analysis B', zorder=3)
ax0.plot(allyrs, mean_rate_ages, color='k', label = 'Analysis D Mean', zorder=3)
ax0.fill_between(allyrs, low_rate2_ages, high_rate2_ages, color = 'lightgray', edgecolor='darkgray', linewidth=0.5, label='Analysis D 1.65$\sigma$', zorder=1)
ax0.fill_between(allyrs, low_rate1_ages, high_rate1_ages, color='darkgray', edgecolor='gray', linewidth=0.5, label = 'Analysis D 1$\sigma$', zorder=1)
ax0.set_xlabel('Years Before Present (YBP)')
ax0.set_ylabel('Storm Frequency ($\lambda$)')
ax0.grid(color='lightgray', linestyle=':', linewidth=0.7, axis='both', which='major', zorder=0)
ax0.set_xlim(-100,8100)
ax0.legend()


# ax1.stem(allyrs, strmyrs, 'k', markerfmt =" ", basefmt=" ")
ax1.plot(allyrs, freq_array_100, label='Analysis A', color='navy',zorder=2)
ax1.plot(allyrs, l50, color='tab:blue', linewidth=2, label='Analysis B', zorder=3)
ax1.plot(allyrs, mean_rate_noages, color='k', label = 'Analysis C Mean', zorder=3)
ax1.fill_between(allyrs, low_rate2_noages, high_rate2_noages, color = 'lightgray', edgecolor='darkgray', linewidth=0.5, label='Analysis C 1.65$\sigma$', zorder=1)
ax1.fill_between(allyrs, low_rate1_noages, high_rate1_noages, color='darkgray', edgecolor='gray', linewidth=0.5, label = 'Analysis C 1$\sigma$', zorder=1)
ax1.set_xlabel('Years Before Present (YBP)')
ax1.set_ylabel('Storm Frequency ($\lambda$)')
ax1.grid(color='lightgray', linestyle=':', linewidth=0.7, axis='both', which='major', zorder=0)
ax1.set_xlim(-100,8100)
ax1.legend()


#%% Calculate Percentage Outside Confidence Intervals

#### Analysis B outside 1 sigma
l50_below_1sig = []
for i in range(0, len(l50)):
    if i >= 3000 and l50[i] < low_rate1_ages[i]:
        l50_below_1sig.append(1)
        
l50_above_1sig = []
for i in range(0, len(freq_array_100)):
    if i >= 3000 and l50[i] > high_rate1_ages[i]:
        l50_above_1sig.append(1)

l50_per_above = np.sum(l50_above_1sig)/len(l50[3000:])*100
l50_per_below = np.sum(l50_below_1sig)/len(l50[3000:])*100

#### Analysis A outside 1 sigma
unif_below_1sig = []
for i in range(0, len(freq_array_100)):
    if freq_array_100[i] < low_rate1_ages[i] and low_rate1_ages[i] > 0.0001:
        unif_below_1sig.append(1)
        
unif_above_1sig = []
for i in range(0, len(freq_array_100)):
    if freq_array_100[i] > high_rate1_ages[i]:
        unif_above_1sig.append(1)

unif_per_above = np.sum(unif_above_1sig)/len(freq_array_100)*100
unif_per_below = np.sum(unif_below_1sig)/len(freq_array_100)*100

#%% Save Bootstrapping Data

np.savetxt(r"JPC05/For_Publication/core_data/l_matrix.csv", l_matrix_noages, delimiter=',')
np.savetxt(r"JPC05/For_Publication/core_data/l_matrix_ages.csv", l_matrix_ages, delimiter=',')
np.savetxt(r"JPC05/For_Publication/core_data/lp_matrix.csv", lp_matrix_noages, delimiter=',')
np.savetxt(r"JPC05/For_Publication/core_data/lp_matrix_ages.csv", lp_matrix_ages, delimiter=',')
# np.savetxt(r"JPC05/core_data/event_matrix_bool.csv", event_matrix_bool, delimiter=',')

#%% Load in Simulation Data if not running now

RMSE_unif_boot = np.genfromtxt(r"JPC05/For_Publication/core_data/RMSE_unif.csv", delimiter=',', skip_header=0)
RMSE_gauss_boot = np.genfromtxt(r"JPC05/For_Publication/core_data/RMSE_gauss.csv",  delimiter=',', skip_header=0)
RMSE_noages_boot = np.genfromtxt(r"JPC05/For_Publication/core_data/RMSE_noages.csv", delimiter=',', skip_header=0)
RMSE_ages_boot = np.genfromtxt(r"JPC05/For_Publication/core_data/RMSE_ages.csv", delimiter=',', skip_header=0)

RMSE_unif2_boot = np.genfromtxt(r"JPC05/For_Publication/core_data/RMSE_unif2.csv", delimiter=',', skip_header=0)
RMSE_gauss2_boot = np.genfromtxt(r"JPC05/For_Publication/core_data/RMSE_gauss2.csv", delimiter=',', skip_header=0)
RMSE_noages2_boot = np.genfromtxt(r"JPC05/For_Publication/core_data/RMSE_noages2.csv", delimiter=',', skip_header=0)
RMSE_ages2_boot = np.genfromtxt(r"JPC05/For_Publication/core_data/RMSE_ages2.csv", delimiter=',', skip_header=0)

RMSE_unif3_boot = np.genfromtxt(r"JPC05/For_Publication/core_data/RMSE_unif3.csv", delimiter=',', skip_header=0)
RMSE_gauss3_boot = np.genfromtxt(r"JPC05/For_Publication/core_data/RMSE_gauss3.csv", delimiter=',', skip_header=0)
RMSE_noages3_boot = np.genfromtxt(r"JPC05/For_Publication/core_data/RMSE_noages3.csv", delimiter=',', skip_header=0)
RMSE_ages3_boot = np.genfromtxt(r"JPC05/For_Publication/core_data/RMSE_ages3.csv", delimiter=',', skip_header=0)


#%% Simulation of Nonstationary Poisson Process

##############################################################################

#Perform the simulation experiment to test analysis performance here
#THIS TAKES 3 DAYS TO RUN
#comment from here to the end of the code if you have already run it

##############################################################################

#%% Simulate Underlying Storm Process

t = np.arange(0,3000) #time of record
lmean = rmean/100 #use average lambda from this record (rmean is in storms per century so divide by 100 to get it in years)
lstd = rstd/100 #convert storms/century to storms/year
tau = 500 #use 500 as the lambda time variation which roughly corresponds to the size of the active interval w/o age uncertainty

#### Simulate Idealized Storm Process
lsim = lmean+lstd*np.sin(2*np.pi*t/tau)


#### Randomly Sample the Underlying Process for storm years
strmyrs_sim = []
for i in range(0,len(lsim)):
    strmyrs_sim.append(np.random.poisson(lsim[i])) #draw from a poisson distribution of the sinusoid
strmyrs_sim_bool=np.array(strmyrs_sim)    #convert to array not list

storm_years = np.argwhere(strmyrs_sim_bool > 0)

#%% Simulate Gaussian Psuedostorm Record

h=50
fake_rate = rate_parameter(t,storm_years,h)
fake_unif = uniform_kernel(t,strmyrs_sim_bool,h)

#Plot Simulated Gaussian Record
plt.figure(dpi=150)
plt.stem(t,strmyrs_sim, 'k', markerfmt =" ", basefmt=" ", label='true storm years')
plt.plot(t,lsim*100,'tab:blue', linewidth=2, label='idealized storm process')
plt.plot(t, fake_rate, 'dimgray', linewidth=2, label='true rate parameter')
# plt.plot(t, fake_unif, 'lightgray', linewidth=2, label='')
plt.grid(color='lightgray', linestyle=':', linewidth=0.7, axis='both', which='major')
plt.ylabel('Storm Frequency ($\lambda$)')
plt.xlabel('YBP')
plt.legend()

#%% Perturb Storm Years

#use average standard deviation of storm deposits in JPCC record
dstd = np.mean(np.std(age_dist,axis=0)) #about 56 years

#pull new storm years from random normal distributions assuming the mean is the simulated storm year and the std is JPCC average event std
strmyrs_pert = []
for i in range(0, len(storm_years)):
    strmyrs_pert.append(int(np.random.normal(storm_years[i],dstd))) #save as integers
storm_years_pert=np.array(strmyrs_pert)

#make age distributions for boostrapping
strmyrs_pert_dist = np.nan*np.zeros((1000, len(storm_years)))
for i in range(0,len(storm_years)):
   strmyrs_pert_dist[:,i] = np.random.normal(storm_years[i],dstd,1000) #use normal distribution even though we know actual distributions aren't normal

#make boolean storm index
strmyrs_pert_bool = np.zeros(len(t))
for i in range(0, len(storm_years_pert)):
    for ii in range(0, len(t)):
        if storm_years_pert[i] == t[ii]:
           strmyrs_pert_bool[ii] = 1 
    

#### plot the disturbed storm years (sanity check)
plt.figure(dpi=150)
plt.stem(t,strmyrs_sim, 'k', markerfmt =" ", basefmt=" ", label='simulated storm years')
plt.stem(t, strmyrs_pert_bool, 'tab:red', markerfmt =" ", basefmt=" ", label='perturbed storm years')
plt.plot(t,lsim*100, 'tab:blue', linewidth=2, label='underlying storm process')
plt.plot(t, fake_rate, 'dimgray', linewidth=2, label='simulated rate parameter')
plt.grid(color='lightgray', linestyle=':', linewidth=0.7, axis='both', which='major')
plt.ylabel('Storm Frequency ($\lambda$)')
plt.xlabel('YBP')
plt.legend()


#%% Reconstruct Psuedostorm Record (median ages)

#### Uniform 
h=50
unif_pert = uniform_kernel(t,strmyrs_pert_bool,h)

#### Gaussian 
gauss_pert = rate_parameter(t,storm_years_pert,h)



#%% Reconstruct Psuedostrom Record with Bootstrapping (without ages)
start_time = time.time()

#set number of iterations
iterations = 1000

#initialize empty matrix
lp_matrix_noages = np.zeros((iterations, len(t)))

for ii in range(0,iterations):
  
    #create a new set of event years with replacement from event years dataset
    rmax = len(storm_years_pert) #maximum index
    rmin = 0 #minimum index
    event_ind = np.random.randint(rmin, high=rmax, size=len(storm_years_pert)) #rmax is exclusive of the final element 
    eventyrs_boot = storm_years_pert[event_ind]
    
    #perform rate parameter analysis on the new timeseries of event years
    l_temp = rate_parameter(t, eventyrs_boot, h)
    
    lp_matrix_noages[ii,:] = l_temp
    
print('time to run bootstrapping is ', (time.time() - start_time)/60, ' minutes')   


#%% Plot Psuedostrom Record with Bootstrapping (without ages)

####Get statistics
mean_prate_noages = np.mean(lp_matrix_noages, axis=0)
low_prate2_noages = np.percentile(lp_matrix_noages, 5, axis=0)
low_prate1_noages = np.percentile(lp_matrix_noages, 16, axis=0)
high_prate2_noages = np.percentile(lp_matrix_noages, 95, axis=0)
high_prate1_noages = np.percentile(lp_matrix_noages, 84, axis=0)

#### Plot
fig,ax = plt.subplots(figsize=(12,4), dpi=150)
# ax.stem(strmyrs_sim, 'k', markerfmt =" ", basefmt=" ", label='simulated storm years') 
ax.plot(lsim*100,'tab:blue', linewidth=2, label='underlying storm process',zorder=2)
# ax.stem(strmyrs_pert_bool, 'b', markerfmt =" ", basefmt=" ", label='perturbed storm years')
ax.plot(t, mean_prate_noages, color='k', label = 'boostrapped no ages',zorder=3)
ax.fill_between(t, low_prate2_noages, high_prate2_noages, color = 'lightgray', edgecolor='darkgray', linewidth=0.5, label='1.65$\sigma$', zorder=1)
ax.fill_between(t, low_prate1_noages, high_prate1_noages, color='darkgray', edgecolor='gray', linewidth=0.5, label = '1$\sigma$', zorder=1)
ax.plot(t, gauss_pert, color='gray', linewidth=1, label='gaussian',zorder=2)
ax.plot(t, unif_pert, color='dimgray', linewidth=1, label='uniform',zorder=2)
ax.grid(color='lightgray', linestyle=':', linewidth=0.7, axis='both', which='major',zorder=0)
ax.legend()
ax.set_ylabel('storm frequency ($\lambda$')
ax.set_xlabel('YBP')

#%% Reconstruct Psuedostrom Record with Bootstrapping (with ages)
start_time = time.time()

#set number of iterations
iterations = 1000

lp_matrix_ages = np.zeros((iterations, len(t)))

for ii in range(0,iterations):
    
    #draw a random event age from the fake probability distributions
    amax = len(strmyrs_pert_dist) #maximum index
    amin = 0 #minimum index
    age_ind = np.random.randint(amin, high=amax, size=1) #get random index
    
    #make alternative set of event years from the simulated distributions
    eventyrs_temp = strmyrs_pert_dist[age_ind,:]
    eventyrs_int = np.squeeze(eventyrs_temp.astype(int)) #convert to integer
    
    #create a new set of event years with replacement
    rmax = len(storm_years_pert) #maximum index
    rmin = 0 #minimum index
    event_ind = np.random.randint(rmin, high=rmax, size=len(storm_years_pert)) #rmax is exclusive of the final element so add 1
    eventyrs_boot = eventyrs_int[event_ind]
    
    #perform rate parameter analysis on the new timeseries of event years
    l_temp = rate_parameter(t, eventyrs_boot, h)
    
    lp_matrix_ages[ii,:] = l_temp


print('time to run bootstrapping is ', (time.time() - start_time)/60, ' minutes')    

#%% Plot Reconstructed Psuedostrom Record with Bootstrapping (with ages)

#### Get statistics
mean_prate_ages = np.mean(lp_matrix_ages, axis=0)
low_prate2_ages = np.percentile(lp_matrix_ages, 5, axis=0)
low_prate1_ages = np.percentile(lp_matrix_ages, 16, axis=0)
high_prate2_ages = np.percentile(lp_matrix_ages, 95, axis=0)
high_prate1_ages = np.percentile(lp_matrix_ages, 84, axis=0)


#### Plot
fig,ax = plt.subplots(figsize=(12,4), dpi=150)
# ax.stem(strmyrs_sim, 'k', markerfmt =" ", basefmt=" ", label='simulated storm years')
ax.plot(lsim*100,'tab:blue', linewidth=2, label='underlying storm process',zorder=2)
# ax.stem(strmyrs_pert_bool, 'b', markerfmt =" ", basefmt=" ", label='perturbed storm years')
ax.plot(t, mean_prate_ages, color='k', linewidth=2, label = 'boostrapped with ages',zorder=3)
ax.fill_between(t, low_prate2_ages, high_prate2_ages, color = 'lightgray', edgecolor='darkgray', linewidth=0.5, label='1.65$\sigma$', zorder=1)
ax.fill_between(t, low_prate1_ages, high_prate1_ages, color='darkgray', edgecolor='gray', linewidth=0.5, label = '1$\sigma$', zorder=1)
ax.plot(t, gauss_pert, color='gray', linewidth=1, label='gaussian',zorder=2)
ax.plot(t, unif_pert, color='dimgray', linewidth=1, label='uniform',zorder=2)
ax.grid(color='lightgray', linestyle=':', linewidth=0.7, axis='both', which='major',zorder=0)
ax.legend()
ax.set_ylabel('storm frequency ($\lambda$')
ax.set_xlabel('YBP')

#%% Get RMSE Between Storm Processes

def rmse(observed, predicted):
    residuals = predicted - observed
    residuals_sq = residuals**2
    mse = np.mean(residuals_sq)
    rmse = np.sqrt(mse)
    
    return(rmse)

#### Underlying Storm Process
#uniform
RMSE_unif = rmse(unif_pert,lsim)
#gaussian
RMSE_gauss = rmse(gauss_pert, lsim)
#bootstrap no ages
RMSE_noages = rmse(mean_prate_noages, lsim)
#bootstrap ages
RMSE_ages = rmse(mean_prate_ages, lsim)

#### Simulated Storm Record
#uniform
RMSE_unif2 = rmse(unif_pert,fake_rate)
#gaussian
RMSE_gauss2 = rmse(gauss_pert, fake_rate)
#bootstrap no ages
RMSE_noages2 = rmse(mean_prate_noages, fake_rate)
#bootstrap ages
RMSE_ages2 = rmse(mean_prate_ages, fake_rate)


#### Reconstructed Storm Record
#uniform
RMSE_unif3 = rmse(unif_pert,fake_unif)
#gaussian
RMSE_gauss3 = rmse(gauss_pert, fake_unif)
#bootstrap no ages
RMSE_noages3 = rmse(mean_prate_noages, fake_unif)
#bootstrap ages
RMSE_ages3 = rmse(mean_prate_ages, fake_unif)

#%% Plot Psudeostorm Simulation with Reconstructions

fig,ax = plt.subplots(4,1, dpi=130, constrained_layout=True)
ax0 = ax[0]
ax1 = ax[1]
ax2 = ax[2]
ax3 = ax[3]

#Plot the Underlying Storm Process
ax0.stem(t,strmyrs_sim, 'k', markerfmt =" ", basefmt=" ", label='simulated storm years')
ax0.plot(t,lsim*100, 'tab:blue', linewidth=2, label='underlying storm process')
ax0.grid(color='lightgray', linestyle=':', linewidth=0.7, axis='both', which='major')
ax0.set_ylabel('Storm Frequency ($\lambda$)')
ax0.set_xlabel('YBP')
ax0.legend()

#Plot the Simulated Psuedostorm Record (Gaussian)
ax1.stem(t,strmyrs_sim, 'k', markerfmt =" ", basefmt=" ", label='simulated storm years')
ax1.stem(t, strmyrs_pert_bool, 'tab:orange', markerfmt =" ", basefmt=" ", label='"reconstructed" storm years')
ax1.plot(t, fake_rate, color='dimgray', linewidth=2, label='simulated frequency (Analysis B)')
ax1.plot(t, gauss_pert, color='tab:red', linewidth=1.5, label='reconstructed frequency (Analysis B)',zorder=2)
ax1.grid(color='lightgray', linestyle=':', linewidth=0.7, axis='both', which='major',zorder=0)
ax1.legend()
ax1.set_ylabel('storm frequency ($\lambda$)')
ax1.set_xlabel('YBP')

#Plot the Simulated Psuedostorm Record (Uniform)
ax2.stem(t,strmyrs_sim, 'k', markerfmt =" ", basefmt=" ", label='simulated storm years')
ax2.stem(t, strmyrs_pert_bool, 'tab:orange', markerfmt =" ", basefmt=" ", label='"reconstructed" storm years')
ax2.plot(t, fake_unif, color='dimgray', linewidth=2, label='simulated frequency (Analysis A)')
ax2.plot(t, unif_pert, color='tab:red', linewidth=1.5, label='reconstructed frequency (Analysis A)',zorder=2)
ax2.grid(color='lightgray', linestyle=':', linewidth=0.7, axis='both', which='major',zorder=0)
ax2.legend()
ax2.set_ylabel('storm frequency ($\lambda$)')
ax2.set_xlabel('YBP')

#Plot the Reconstructed Storm Frequency
ax3.plot(lsim*100,'tab:blue', linewidth=2, label='underlying storm process',zorder=2)
ax3.plot(t, fake_rate, 'darkblue', linewidth=2, label='simulated frequency (Analysis B)')
ax3.plot(t, mean_prate_ages, color='k', linewidth=2, label = 'Reconstructed Frequency (Analysis D Mean)',zorder=3)
ax3.fill_between(t, low_prate2_ages, high_prate2_ages, color = 'lightgray', edgecolor='darkgray', linewidth=0.5, label='Analysis D 1.65$\sigma$', zorder=1)
ax3.fill_between(t, low_prate1_ages, high_prate1_ages, color='darkgray', edgecolor='gray', linewidth=0.5, label = 'Analysis D 1$\sigma$', zorder=1)
# ax2.plot(t, gauss_pert, color='gray', linewidth=1, label='Analysis B',zorder=2)
# ax2.plot(t, unif_pert, color='dimgray', linewidth=1, label='Analysis A',zorder=2)
ax3.grid(color='lightgray', linestyle=':', linewidth=0.7, axis='both', which='major',zorder=0)
ax3.legend()
ax3.set_ylabel('storm frequency ($\lambda$)')
ax3.set_xlabel('YBP')



#%% Bootstrap Psuedostorm Simulation

start_time = time.time()

# outer loop is the number of simulations (only do 100 simulations but bootstrap 1000 times within the sim)
iterations1 = 100

#### initialize empty arrays
RMSE_unif_boot = np.nan * np.zeros(iterations1)
RMSE_gauss_boot = np.nan * np.zeros(iterations1)
RMSE_noages_boot = np.nan * np.zeros(iterations1)
RMSE_ages_boot = np.nan * np.zeros(iterations1)

RMSE_unif2_boot = np.nan * np.zeros(iterations1)
RMSE_gauss2_boot = np.nan * np.zeros(iterations1)
RMSE_noages2_boot = np.nan * np.zeros(iterations1)
RMSE_ages2_boot = np.nan * np.zeros(iterations1)

RMSE_unif3_boot = np.nan * np.zeros(iterations1)
RMSE_gauss3_boot = np.nan * np.zeros(iterations1)
RMSE_noages3_boot = np.nan * np.zeros(iterations1)
RMSE_ages3_boot = np.nan * np.zeros(iterations1)

faker_Gauss_outside_ages = np.nan * np.zeros(iterations1)
faker_Unif_outside_ages = np.nan * np.zeros(iterations1)
lsim_outside_ages = np.nan * np.zeros(iterations1)

faker_Gauss_outside_noages = np.nan * np.zeros(iterations1)
faker_Unif_outside_noages = np.nan * np.zeros(iterations1)
lsim_outside_noages = np.nan * np.zeros(iterations1)

iterations = 1000

#### Monte Carlo Simulation
for nn in range(0, iterations1):

    start_time1 = time.time()
    
    #randomly choose storm years for the record
    strmyrs_sim = []
    for i in range(0,len(lsim)):
        strmyrs_sim.append(np.random.poisson(lsim[i])) #draw from a the underlying storm process
    strmyrs_sim_bool=np.array(strmyrs_sim)    #convert to array not list
    
    #get storm years
    storm_years = np.argwhere(strmyrs_sim_bool > 0)

    #simulate psuedostorm frequency
    h=50
    fake_rate = rate_parameter(t,storm_years,h) #Gaussian
    fake_unif = uniform_kernel(t, strmyrs_sim_bool, h) #uniform

    #perturb storm years by sampling from random normal distributions assuming the mean is the simulated storm year and the std is this record std
    strmyrs_pert = []
    for i in range(0, len(storm_years)):
        strmyrs_pert.append(int(np.random.normal(storm_years[i],dstd))) #save as integers
    storm_years_pert=np.array(strmyrs_pert)

    #make age distributions for boostrapping
    strmyrs_pert_dist = np.nan*np.zeros((1000, len(storm_years)))
    for i in range(0,len(storm_years)):
       strmyrs_pert_dist[:,i] = np.random.normal(storm_years[i],dstd,1000)

    #make boolean storm index
    strmyrs_pert_bool = np.zeros(len(t))
    for i in range(0, len(storm_years_pert)):
        for ii in range(0, len(t)):
            if storm_years_pert[i] == t[ii]:
               strmyrs_pert_bool[ii] = 1 
               
    #reconstruct storm frequency 
    h=50
    unif_pert = uniform_kernel(t,strmyrs_pert_bool,h) #uniform
    gauss_pert = rate_parameter(t,storm_years_pert,h) #gauss
    
    
#### Reconstruct Psuedostorm Frequency with Bootstrapping (without ages)
    #initialize empty matrix
    lp_matrix_noages = np.zeros((iterations, len(t)))
    
    # start_b1 = time.time()
    for ii in range(0,iterations):
        
      
        #create a new set of event years with replacement from event years dataset
        rmax = len(storm_years_pert) #maximum index
        rmin = 0 #minimum index
        event_ind = np.random.randint(rmin, high=rmax, size=len(storm_years_pert)) #rmax is exclusive of the final element 
        eventyrs_boot = storm_years_pert[event_ind]
        
        #perform rate parameter analysis on the new timeseries of event years
        l_temp = rate_parameter(t, eventyrs_boot, h)
        
        lp_matrix_noages[ii,:] = l_temp
        
        #Get statistics
        mean_prate_noages = np.mean(lp_matrix_noages, axis=0)
        low_prate2_noages = np.percentile(lp_matrix_noages, 5, axis=0)
        low_prate1_noages = np.percentile(lp_matrix_noages, 16, axis=0)
        high_prate2_noages = np.percentile(lp_matrix_noages, 95, axis=0)
        high_prate1_noages = np.percentile(lp_matrix_noages, 84, axis=0)
    
    # print('time to run no ages bootstrapping is ', (time.time() - start_b1)/60, ' minutes')      
    
#### Reconstruct Psuedostorm Frequency with Bootstrapping (with ages)
    #initialize empty matrix
    lp_matrix_ages = np.zeros((iterations, len(t)))
    
    # start_b2 = time.time()
    for ii in range(0,iterations):
        
        #draw a random event age from the Bacon probability distributions
        amax = len(strmyrs_pert_dist) #maximum index
        amin = 0 #minimum index
        age_ind = np.random.randint(amin, high=amax, size=1) #get random index
        
        #make alternative set of event years from the random age model iteration
        eventyrs_temp = strmyrs_pert_dist[age_ind,:]
        eventyrs_int = np.squeeze(eventyrs_temp.astype(int)) #convert to integer
        
        #create a new set of event years with replacement
        rmax = len(storm_years_pert) #maximum index
        rmin = 0 #minimum index
        event_ind = np.random.randint(rmin, high=rmax, size=len(storm_years_pert)) #rmax is exclusive of the final element so add 1
        eventyrs_boot = eventyrs_int[event_ind]
        
        #perform rate parameter analysis on the new timeseries of event years
        l_temp = rate_parameter(t, eventyrs_boot, h)
        
        lp_matrix_ages[ii,:] = l_temp
        
        #Get statistics
        mean_prate_ages = np.mean(lp_matrix_ages, axis=0)
        low_prate2_ages = np.percentile(lp_matrix_ages, 5, axis=0)
        low_prate1_ages = np.percentile(lp_matrix_ages, 16, axis=0)
        high_prate2_ages = np.percentile(lp_matrix_ages, 95, axis=0)
        high_prate1_ages = np.percentile(lp_matrix_ages, 84, axis=0)
        
    # print('time to run ages bootstrapping is ', (time.time() - start_b2)/60, ' minutes')   
    
#### GET RMSE

    ### Underlying Storm Process
    RMSE_unif_boot[nn] = rmse(unif_pert,lsim*100)
    #gaussian
    RMSE_gauss_boot[nn] = rmse(gauss_pert, lsim*100)
    #bootstrap no ages
    RMSE_noages_boot[nn] = rmse(mean_prate_noages, lsim*100)
    #bootstrap ages
    RMSE_ages_boot[nn] = rmse(mean_prate_ages, lsim*100)
    
    ### Simulated Gaussian Frequency
    #uniform
    RMSE_unif2_boot[nn] = rmse(unif_pert,fake_rate)
    #gaussian
    RMSE_gauss2_boot[nn] = rmse(gauss_pert, fake_rate)
    #bootstrap no ages
    RMSE_noages2_boot[nn] = rmse(mean_prate_noages, fake_rate)
    #bootstrap ages
    RMSE_ages2_boot[nn] = rmse(mean_prate_ages, fake_rate)
    
    ### Simulated Uniform Frequency
    #uniform
    RMSE_unif3_boot[nn] = rmse(unif_pert,fake_unif)
    #gaussian
    RMSE_gauss3_boot[nn] = rmse(gauss_pert, fake_unif)
    #bootstrap no ages
    RMSE_noages3_boot[nn] = rmse(mean_prate_noages, fake_unif)
    #bootstrap ages
    RMSE_ages3_boot[nn] = rmse(mean_prate_ages, fake_unif)
    
    
  #### Estimate how often analysis falls outside confidence intervals from Analysis D
    lsim_out = []
    faker_gauss_out = []
    faker_unif_out = []
    
    ###underlying storm process
    for ii in range(0, len(lsim)):
        if lsim[ii]*100 < low_prate2_ages[ii] or lsim[ii]*100 > high_prate2_ages[ii]:
            lsim_out.append(1)
    
    lsim_outside_ages[nn] = np.sum(lsim_out)/len(lsim)*100
    
    ###simulated gaussian psuedostorm frequency
    for ii in range(0,len(fake_rate)):
        if fake_rate[ii] < low_prate2_ages[ii] or fake_rate[ii] > high_prate2_ages[ii]:
            faker_gauss_out.append(1)
    
    faker_Gauss_outside_ages[nn] = np.sum(faker_gauss_out)/len(fake_rate)*100
    
    ###simulated uniform psuedostorm frequency
    for ii in range(0,len(fake_unif)):
        if fake_unif[ii] < low_prate2_ages[ii] or fake_unif[ii] > high_prate2_ages[ii]:
            faker_unif_out.append(1)
    
    faker_Gauss_outside_ages[nn] = np.sum(faker_unif_out)/len(fake_unif)*100
    
    
    
  #### Estimate how often analysis falls outside confidence intervals from Analysis C
    lsim_out = []
    faker_gauss_out = []
    faker_unif_out = []
    
    ###underlying storm process
    for ii in range(0, len(lsim)):
        if lsim[ii]*100 < low_prate2_noages[ii] or lsim[ii]*100 > high_prate2_noages[ii]:
            lsim_out.append(1)
    
    lsim_outside_noages[nn] = np.sum(lsim_out)/len(lsim)*100
    
    ###simulated gaussian psuedostorm frequency
    for ii in range(0,len(fake_rate)):
        if fake_rate[ii] < low_prate2_noages[ii] or fake_rate[ii] > high_prate2_noages[ii]:
            faker_gauss_out.append(1)
    
    faker_Gauss_outside_noages[nn] = np.sum(faker_gauss_out)/len(fake_rate)*100
    
    ###simulated uniform psuedostorm frequency
    for ii in range(0,len(fake_unif)):
        if fake_unif[ii] < low_prate2_noages[ii] or fake_unif[ii] > high_prate2_noages[ii]:
            faker_unif_out.append(1)
    
    faker_Gauss_outside_noages[nn] = np.sum(faker_unif_out)/len(fake_unif)*100
    
        
    print('time to run one bootstrap is ', (time.time() - start_time1)/60, ' minutes')  
    print('now on iteration', nn) 
    
print('time to run bootstrapping is ', (time.time() - start_time)/60, ' minutes')  


    
#%% Save Simulation Data 

np.savetxt(r"JPC05/For_Publication/core_data/RMSE_unif.csv", RMSE_unif_boot, delimiter=',')
np.savetxt(r"JPC05/For_Publication/core_data/RMSE_gauss.csv", RMSE_gauss_boot, delimiter=',')
np.savetxt(r"JPC05/For_Publication/core_data/RMSE_noages.csv", RMSE_noages_boot, delimiter=',')
np.savetxt(r"JPC05/For_Publication/core_data/RMSE_ages.csv", RMSE_ages_boot, delimiter=',')

np.savetxt(r"JPC05/For_Publication/core_data/RMSE_unif2.csv", RMSE_unif2_boot, delimiter=',')
np.savetxt(r"JPC05/For_Publication/core_data/RMSE_gauss2.csv", RMSE_gauss2_boot, delimiter=',')
np.savetxt(r"JPC05/For_Publication/core_data/RMSE_noages2.csv", RMSE_noages2_boot, delimiter=',')
np.savetxt(r"JPC05/For_Publication/core_data/RMSE_ages2.csv", RMSE_ages2_boot, delimiter=',')

np.savetxt(r"JPC05/For_Publication/core_data/RMSE_unif3.csv", RMSE_unif3_boot, delimiter=',')
np.savetxt(r"JPC05/For_Publication/core_data/RMSE_gauss3.csv", RMSE_gauss3_boot, delimiter=',')
np.savetxt(r"JPC05/For_Publication/core_data/RMSE_noages3.csv", RMSE_noages3_boot, delimiter=',')
np.savetxt(r"JPC05/For_Publication/core_data/RMSE_ages3.csv", RMSE_ages3_boot, delimiter=',')

np.savetxt(r"JPC05/For_Publication/core_data/lsim_outside_ages.csv", lsim_outside_ages, delimiter=',')
np.savetxt(r"JPC05/For_Publication/core_data/faker_Gauss_outside_ages.csv", faker_Gauss_outside_ages, delimiter=',')
np.savetxt(r"JPC05/For_Publication/core_data/faker_Unif_outside_ages.csv", faker_Unif_outside_ages, delimiter=',')

np.savetxt(r"JPC05/For_Publication/core_data/lsim_outside_noages.csv", lsim_outside_noages, delimiter=',')
np.savetxt(r"JPC05/For_Publication/core_data/faker_Gauss_outside_noages.csv", faker_Gauss_outside_noages, delimiter=',')
np.savetxt(r"JPC05/For_Publication/core_data/faker_Unif_outside_noages.csv", faker_Unif_outside_noages, delimiter=',')

print('time to run code is ', (time.time() - start_code)/60, ' minutes')    
