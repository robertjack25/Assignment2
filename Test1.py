"""
Created on Thu Oct 20 15:23:50 2022
@author: robert jack
"""

import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt



data = np.loadtxt( 'robert_ecg.dat' )

print (data)
time = data[:,0]
first = data[:,1]
second = data[:,2]   #channel we selected
third =data[:,3]

ecg_dat = data[:,2]
smaller_sample = ecg_dat[0:2001]

fs = 1000 #sampling frequency [Hz]
Ts = 1/fs
time1 = np.zeros(len(time))
for i in range(len(time)):
    time1[i] = Ts * time[i]

#freq res = fs/M

def filterDesign (fs, f1, f2, M):   #this function creates the bandstop impulse response, h
    
    k1 = int(f1/fs * M)
    k2 = int(f2/fs * M)
    X = np.ones(M)
    X[k1:k2+1] = 0
    X[M-k2:M-k1+1] = 0
    x = np.fft.ifft(X)
    x = np.real(x)
    h = np.zeros(M)
    h[0:int(M/2)] = x[int(M/2):M]
    h[int(M/2):M] = x[0:int(M/2)]
    h = h*np.hamming(M)
    return h

##for removing the DC do we just remove X[0]?

remove_50Hz = filterDesign(1000, 45, 55, 1000)    
arr_length =  len(remove_50Hz)
print('It is: %i' %arr_length)

#letter = 'original'
#fs,audio_dat = wavfile.read('%s.wav' %letter)

class FIR_filter:
    def __init__(self,coeff_vals):    #coeff_vals here is the returned impulse response from filter design
        self.num_taps = len(coeff_vals)
        self.coefficients = coeff_vals
       # self.buffer = [0]*self.num_taps #creates scalar array 
        self.buffer = np.zeros(self.num_taps)

    def dofilter(self,ecgsig_val):  
        self.buffer[0] = ecgsig_val   #why does it need the [0]? initially it had [0]
    
          #y(n) = x(n) * h(n)
        filtered_ecgsig = np.zeros(self.num_taps)
        for n in range(self.num_taps): 
             filtered_ecgsig += self.buffer[n]*self.coefficients[n]  #y(n) = x(n) * h(n)
             self.buffer = np.roll(self.buffer,1)
        return filtered_ecgsig
            
        #need to make this realtime and feed data (ecgsig_val), sample by sample into the function   
        #also, why is the length of the buffer not equal to the length of the taps?
         
            
    
    
######REMEMBER TO ADD IN HIGHPASS FILTERING TO REMOVE BASELINE WANDER    
#filt_ecg = np.zeros(len(smaller_sample))
filt_ecg = []
for n in range(len(smaller_sample)):
    using_FIR = FIR_filter(remove_50Hz) #instantiating the class with impulse response to remove 50Hz
    y = using_FIR.dofilter(smaller_sample[n])
    filt_ecg[n] = y
    
    print(n)
    
    
#-------------------------------------
scipy_filtered = signal.lfilter(remove_50Hz,1,ecg_dat)

#plt.plot(time1, second)
#plt.plot(filterDesign() # removed 50Hz

fig1 = plt.figure(1)
plt.plot(time1, ecg_dat)
plt.xlabel('Time (s)')
plt.ylabel('Voltage')
plt.title('Initial ECG data')
plt.grid(which='major', color='black', linewidth=0.8)
plt.grid(which='minor', color='black', linewidth=0.5)
plt.minorticks_on()
plt.show() 

fig2 = plt.figure(2)
plt.plot(remove_50Hz)
#plt.xlabel('Time (s)')
#plt.ylabel('')
plt.title('Impulse Response')
plt.grid(which='major', color='black', linewidth=0.8)
plt.grid(which='minor', color='black', linewidth=0.5)
plt.minorticks_on()
plt.show()
 

fig3 = plt.figure(3)
plt.plot(filt_ecg)
plt.xlabel('Time (s)')
plt.ylabel('Voltage')
plt.title('Filtered ECG')
plt.grid(which='major', color='black', linewidth=0.8)
plt.grid(which='minor', color='black', linewidth=0.5)
plt.minorticks_on()
plt.show()

fig4 = plt.figure(4)
plt.plot(scipy_filtered)
plt.xlabel('Time (s)')
plt.ylabel('Voltage')
plt.title('Scipy Filtered ECG')
plt.grid(which='major', color='black', linewidth=0.8)
plt.grid(which='minor', color='black', linewidth=0.5)
plt.minorticks_on()
plt.show()
