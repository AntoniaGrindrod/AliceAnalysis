#Program to 
import numpy as np
import matplotlib.pyplot as plt
import os
import csv
from scipy.optimize import curve_fit


##Read in o32 to a 3d array (each element is an event; each element in event is an array containing the data points for the 3 channels)
filename = 'all_oscilloscope_2242.o32'
file = open(filename)

data_arr = []

for i in range(27): #change 27 to the number of events
    data_arr.append([])
    for j in range(4):
        file.readline()
    for j in range(1000):
        line = file.readline()
        line.strip('\n')
        data_arr[i].append((line.split(',')))
        for s in range(len(data_arr[i][j])):
            data_arr[i][j][s] = float(data_arr[i][j][s])
        #data_arr[i][j] = np.array(data_arr[i][j])
#for s in range(len(data_arr[0])):
    #print(data_arr[0][s])
data_arr_new=[]
#Need to make this array so that each event rather has 3 entries (channel) (i.e. we're switching the second two dimensions of the array)
for s in range(len(data_arr)):
    data_arr_new.append([[],[],[]])
    for t in range(len(data_arr[s])):#number of data points in each event
        data_arr_new[s][0].append(data_arr[s][t][0]) #channel 1
        data_arr_new[s][1].append(data_arr[s][t][1]) #channel 2
        data_arr_new[s][2].append(data_arr[s][t][2]) #trigger pulse
    
##Get Time Difference

#Function definition to read in scintillator files into a numpy array

##Arrays

peakPosA1 = []
peakPosA2 = []

peakPosB1 = []
peakPosB2 = []

peakPosC1 =[]
peakPosC2=[]

bad_data=0

## Function for finding roots of linear fit done later
def solve_for_y(poly_coeffs, y):
    pc = poly_coeffs.copy()
    pc[-1] -= y
    return np.roots(pc)

## For loop
for r in range(len(data_arr_new)):
    if True:
        pulse = data_arr_new[r]
        
        
        ## Splitting the Pulses
        
        pulse1 = np.array(pulse[0])# channel 1 data
        pulse2 = np.array(pulse[1])# channel 2
        print(pulse1)
        
        ##Time array as time-step of 4 ns
        t=np.arange(0,len(pulse1)*4,4)
        
        ##Isolating signal Peak Areas
        
        reduced_pulse1=pulse1[350:450]
        reduced_pulse2=pulse2[350:450]
        

        ## Reducing Time arrray to help with indexing
        
        reduced_t=t[350:450]
        ## Turning Points of Signals
        
        peak_1 = np.amin(reduced_pulse1)
        index_1 = np.where(reduced_pulse1 == peak_1)
        
        peak_2 = np.amin(reduced_pulse2)
        index_2 = np.where(reduced_pulse2 == peak_2)
        
        ## Finding small fluctuations away from the zero line
        signal_fluctuations_1 = pulse1[0:350]
        signal_fluctuations_2 = pulse2[0:350]
    
        
        for i in range (len(signal_fluctuations_1)):
            if (signal_fluctuations_1[i] != 0) & (signal_fluctuations_1[i] > 0):
                signal_fluctuations_1 = -signal_fluctuations_1[i]
                break
            
            elif (signal_fluctuations_1[i] != 0) & (signal_fluctuations_1[i] < 0):
                signal_fluctuations_1 = signal_fluctuations_1[i]
                break
        
        if np.isscalar(signal_fluctuations_1)==False:
            signal_fluctuations_1=0

        for i in range (len(signal_fluctuations_2)):
            if (signal_fluctuations_2[i] != 0) & (signal_fluctuations_2[i] > 0):
                signal_fluctuations_2 = -signal_fluctuations_2[i]
                break
            
            elif (signal_fluctuations_2[i] != 0) & (signal_fluctuations_2[i] < 0):
                signal_fluctuations_2 = signal_fluctuations_2[i]
                break
        
        if np.isscalar(signal_fluctuations_2)==False:
            signal_fluctuations_2=0
                
        ## Checking to see if there are any "bad" data sets
        pulse_test_1 = pulse1
        pulse_test_2 = pulse2

        for i in range(len(pulse_test_1)):
            if (np.abs(pulse_test_1[i])==np.abs(signal_fluctuations_1)):
                pulse_test_1[i]=0
        
        for i in range(len(pulse_test_2)):
            if (np.abs(pulse_test_2[i])==np.abs(signal_fluctuations_2)):
                pulse_test_2[i]=0
                
        if (np.all(pulse_test_1==0)==True) or (np.all(pulse_test_2==0)==True):
            peakPosA1.append('NaN')
            peakPosA2.append('NaN')
            
            peakPosB1.append('NaN')
            peakPosB2.append('NaN')
            
            peakPosC1.append('NaN')
            peakPosC2.append('NaN')
            
            bad_data=1
        
        else:
            
            ##Method 1: measuring time difference between the minima
            if (len(reduced_t[index_1]) > 1):
                time_1=sum(reduced_t[index_1])/len(reduced_t[index_1])
            else:
                time_1=reduced_t[index_1][0]
            peakPosA1.append(time_1)
        
            if (len(reduced_t[index_2]) >1):
                time_2=sum(reduced_t[index_2])/len(reduced_t[index_2])
            else:
                time_2=reduced_t[index_2][0]
            peakPosA2.append(time_2)
    
            ## Method 2: Threshold value done in 2020
            diff_1 = peak_1/2
            count = 0

            for i in range(2,len(reduced_t)):
                if np.abs(reduced_pulse1[i-2] - reduced_pulse1[i]) > diff_1:
                    count = i-2
                    pos = reduced_t[count]
                    peakPosB1.append(pos)  
                    break
                    
            diff_2 = peak_2/2
            count = 0

            for i in range(2,len(reduced_t)):
                if np.abs(reduced_pulse2[i-2] - reduced_pulse2[i]) > diff_2:
                    count = i-2
                    pos = reduced_t[count]
                    peakPosB2.append(pos)  
                    break
            
            ## Method 3: Dunno what to call it? Linear fit for non-discrete starting value of the pulse
            
            ## Preliminary set-up of arrays
            fit_data1=[]
            fit_t1=[]
            for i in range (len(reduced_pulse1)):
                if ((reduced_pulse1[i] < signal_fluctuations_1) & (reduced_pulse1[i] >= peak_1)):
                    fit_data1.append(reduced_pulse1[i-1])
                    fit_t1.append(reduced_t[i-1])
                    
                    if (reduced_pulse1[i] == peak_1):
                        fit_data1.append(reduced_pulse1[i])
                        fit_t1.append(reduced_t[i])
                        break
            
            fit_data2=[]
            fit_t2=[]
            for i in range (len(reduced_pulse2)):
                if ((reduced_pulse2[i] < signal_fluctuations_2) & (reduced_pulse2[i] >= peak_2)):
                    fit_data2.append(reduced_pulse2[i-1])
                    fit_t2.append(reduced_t[i-1])
                    
                    if (reduced_pulse2[i] == peak_2):
                        fit_data2.append(reduced_pulse2[i])
                        fit_t2.append(reduced_t[i])
                        break
            
            ## Fit of the linear function
            fit_1 = np.polyfit(fit_t1,fit_data1,1)
            fit_2 = np.polyfit(fit_t2,fit_data2,1)
            
            if signal_fluctuations_1 >= signal_fluctuations_2:
                eval = signal_fluctuations_2
            else:
                eval = signal_fluctuations_1
            
            peakPosC1.append(solve_for_y(fit_1,2*eval)[0].astype(float))
            peakPosC2.append(solve_for_y(fit_2,2*eval)[0].astype(float))
            


## Converting to numpy arrays for ease of use
peakPosA1=np.array(peakPosA1)
peakPosA2=np.array(peakPosA2)

peakPosB1=np.array(peakPosB1)
peakPosB2=np.array(peakPosB2)

peakPosC1=np.array(peakPosC1)
pealPosC2=np.array(peakPosC2)

##Removing any points from "bad" data sets
if (bad_data == 1):
    peakPosA1 = np.delete(peakPosA1,np.where(peakPosA1 == 'NaN')).astype(float)
    peakPosA2 = np.delete(peakPosA2,np.where(peakPosA2 == 'NaN')).astype(float)
    peakPosA1=peakPosA1[~np.isnan(peakPosA1)]
    peakPosA2=peakPosA2[~np.isnan(peakPosA2)]

    peakPosB1 = np.delete(peakPosB1,np.where(peakPosB1 == 'NaN')).astype(float)
    peakPosB2 = np.delete(peakPosB2,np.where(peakPosB2 == 'NaN')).astype(float)
    peakPosB1=peakPosB1[~np.isnan(peakPosB1)]
    peakPosB2=peakPosB2[~np.isnan(peakPosB2)]

    peakPosC1 = np.delete(peakPosC1,np.where(peakPosC1 == 'NaN')).astype(float)
    peakPosC2 = np.delete(peakPosC2,np.where(peakPosC2 == 'NaN')).astype(float)
    peakPosC1=peakPosC1[~np.isnan(peakPosC1)]
    peakPosC2=peakPosC2[~np.isnan(peakPosC2)]
    

##Time difference from method 1
tdArr1 = peakPosA1 - peakPosA2

##Time difference from method 2
tdArr2 = peakPosB1 - peakPosB2

##Time difference from method 3
tdArr3 = peakPosC1 - peakPosC2






##Find Angle
distances_arr = [2.865]
distances = np.array(distances_arr) #array to store distances between scintillators (for use in trial to determine which is optimal distance)
u_distance = np.sqrt(2)*(0.001)/(2*np.sqrt(6)) #set this value to the uncertainty on the distance was before: u_distance = (0.001)/(2*np.sqrt(6))
time_res = 5e-9 #use 5ns for trials - maybe change later to actual time res values
u_time_diff = np.sqrt(2)*time_res #same uncertainty on measured time difference in each case, since this is just propagation from time resolution
muon_speed = 0.994*3e8
event_count= np.array([]) #Number of coincident events detected in each trial

#Create 2-d array of scintillator data (each entry is an array storing the scintillator data for each trial)
#put scintillator 1 and 2 data into this array (entry in each corresponsing position must correspond to same event)
time_differences = [] #This must be the scintillator on the "left (further from TRD)" #This must be the scintillator on the "right (closer to TRD)"
##These arrays must store floats or integers, or code won't calculate the difference in their entries properly



#Create 2-d arrays (an array of d arrays, where d is the number of entries in distances array to store theta data for each trial -DON'T ACTUALLY NEED 2d ARRAYS, but had adapted code from distances_trial
theta_data = []
u_theta_data = []

count_failed = [] #stores number of data points meeting condition in line 57
theta_zero_count =[]
count_uncertainty_failed = []
for i in range(len(distances)):
    theta_data.append([])
    u_theta_data.append([])
    time_differences.append(tdArr3) #this is hard-coded, since we only have one distance
    count_failed.append(0)
    theta_zero_count.append(0)
    count_uncertainty_failed.append(0)
    
    
#initialise variables
vt = 0
u_vt = 0
cos_theta = 0
theta = 0

for j in range(len(distances)):
    distance = distances[j]
    print("Distance: ",distance)
    for i in range(len(time_differences[j])):
        time_diff =  np.abs(time_differences[j][i])*1e-9 #differences in time of event as measured by scintillators (array has values in nanoseconds and we want in second)
        vt = muon_speed*time_diff
        
        u_vt = muon_speed*u_time_diff
        cos_theta = (vt)/distance
        if cos_theta>1: #Removing data points that don't fit the model (this might be because they aren't actually muon shower events). In any case, WE CAN'T GET AN ANGLE FROM THESE EVENTS
            count_failed[j] = count_failed[j]+1
        else:
            theta = np.arccos(cos_theta)
            if vt==0: #separate this case since otherwise get an infinite uncertainty
                theta_zero_count[j] = theta_zero_count[j]+1       
            else:
                u_cos_theta = np.sqrt((u_vt/vt)**2+(u_distance/distance)**2)
                u_theta = u_cos_theta/np.sin(theta)
                
                if time_differences[j][i]<0: #If time of event as measured by 1 is before that of 2, the angle is measured clockwise from the horizontal (so our angle will be whatever we get with the calculation subtracted from 180 deg) otherwise it is measured counterclockwise
                    theta = 180 - theta #note the uncertainty on theta will remain the same
                
                if u_theta>5: #Remove angles with uncertainties which are clearly outliers (this condition is rather arbitrary and should be changed - i.e the reason why these data points are giving very high uncertainties should be explored, but we did not have time to do this)
                    count_uncertainty_failed[j] = count_uncertainty_failed[j]+1
                else:
                    theta_data[j].append(theta)
                    u_theta_data[j].append(u_theta)
    
    theta_data[j] = np.array(theta_data[j])
    u_theta_data[j] = np.array(u_theta_data[j])  
        


u_theta_avg = []
for i in range(len(theta_data)):
    u_theta_avg.append(np.mean(u_theta_data[i]))

print("Separation Distance:\t Average Angular Resolution:")
for i in range(len(distances)):
    print(distances[i],'\t',np.mean(u_theta_data[i]))
    
#for i in range(len(distances)):
    #print(distances[i])
    #print(theta_data[i])
    #print(u_theta_data[i])
    #print(count_failed[i])
    #print(count_uncertainty_failed[i]) 
        
