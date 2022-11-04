
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Oct 20 15:23:50 2022
@author: robert jack
"""

import numpy as np
import matplotlib.pyplot as plt


file = 'ecg1.dat'

def CreateImpulseResponse (fs, f1, f2, M):   #this function creates the bandstop impulse response, h
    
    f_dc = 0.5     #remove dc contributions below 0.5Hz
    f_high = 100   #this removes noise above possible ecg frequencies
    bin_dc = int((f_dc/fs) * M)
    bin_f1 = int((f1/fs) * M)
    bin_f2 = int((f2/fs) * M)
    bin_high = int((f_high/fs) * M)
    X = np.ones(M)
    X[:bin_dc] = 0
    X[-bin_dc:M] = 0
    X[bin_f1:bin_f2+1] = 0
    X[-bin_f2:-bin_f1+1] = 0
    X[bin_high:int(M/2)+1] = 0    #made M/2 explicitly an int as error was created otherwise
    X[int(M/2):-bin_high+1] = 0  
    x = np.fft.ifft(X)
    x = np.real(x)
    h = np.zeros(M)
    h[:int(M/2)] = x[int(M/2):M]
    h[int(M/2):M] = x[:int(M/2)]
    h = h*np.hamming(M)
    return h, X

class FIRfilter:
    
    def __init__(self,coeff_vals):    #coeff_vals here is array of returned impulse response from filter design
        self.num_taps = len(coeff_vals)
        self.coefficients = coeff_vals
        self.buffer = np.zeros(self.num_taps)
        self.delta_h = np.zeros(self.num_taps)

    def dofilter(self,ecgsig_val):
        self.buffer = np.roll(self.buffer,1)
        self.buffer[0] = ecgsig_val
        result = np.inner(self.buffer, self.coefficients)
        return result
            
    def doFilterAdaptive(self,signal,noise,learningRate):  #signal is the ecg data(the noisy signal)
        remover = self.dofilter(noise)
        output_e = signal - remover
        self.delta_h = output_e * learningRate * self.buffer
        self.coefficients = np.add(self.coefficients,self.delta_h)
        return output_e


data = np.loadtxt( '%s' %file)
time = data[:,0]
first = data[:,1]
second = data[:,2]   #channel we selected
third = data[:,3]

ecg_dat = data[:,2]
ecg_dat = 1.325*(ecg_dat - 2**23)/2**23 #scaling ADC bit value
#multiplying by 1.325 as this is the magnitude of the smallest and largest possible voltage signal

fs = 1000 #sampling frequency [Hz]
Ts = 1/fs #sampling period [s]

faxis = np.linspace(0,fs,len(ecg_dat))
taxis = np.linspace(0,len(ecg_dat)/fs, len(ecg_dat))

#freq res = fs/M         
fnoise = 50
fnoise2 = 210

impulse_response, X = CreateImpulseResponse(1000, 45, 55, 2000)    

filt_ecg = np.zeros(len(ecg_dat))

using_FIR = FIRfilter(impulse_response)  #instatiating the class

inp_noise = np.zeros(len(ecg_dat))
for n in range(len(ecg_dat)):
    filt_ecg[n] = using_FIR.dofilter(ecg_dat[n]) #using filter to remove 50Hz, DC amd high frequencies

plt.figure()
plt.plot(taxis, ecg_dat)
plt.xlabel('Time (s)')
plt.ylabel('Voltage (V)')
plt.title('Initial ECG data')
plt.grid(which='major', color='black', linewidth=0.8)
plt.grid(which='minor', color='black', linewidth=0.5)
plt.minorticks_on()
plt.savefig('initial_ecg_%s.svg' %file, format = 'svg')

plt.figure()
plt.plot(taxis, filt_ecg)
plt.xlabel('Time(s)')
plt.ylabel('Voltage (V)')
plt.title('Filtered Sitting ECG data')
plt.grid(which='major', color='black', linewidth=0.8)
plt.grid(which='minor', color='black', linewidth=0.5)
plt.minorticks_on()
plt.savefig('ecg_filtered_%s.svg' %file, format = 'svg')

plt.figure()
plt.plot(taxis, ecg_dat)
plt.xlabel('Time(s)')
plt.xlim(8.7,9.5)
plt.ylabel('Voltage (V)')
plt.title('Zoomed In Filtered Sitting ECG data')
plt.grid(which='major', color='black', linewidth=0.8)
plt.grid(which='minor', color='black', linewidth=0.5)
plt.minorticks_on()
plt.savefig('ecg_filtered_ZoomedIn_%s.svg' %file, format = 'svg')

plt.figure()
plt.plot(taxis, filt_ecg)
plt.xlabel('Time(s)')
plt.xlim(9,9.9)
plt.ylabel('Voltage (V)')
plt.title('Zoomed In Filtered ECG data')
plt.grid(which='major', color='black', linewidth=0.8)
plt.grid(which='minor', color='black', linewidth=0.5)
plt.minorticks_on()
plt.savefig('Zoomed_ecg_filtered_%s.svg' %file, format = 'svg')

plt.figure()
#fig6 = plt.figure(6)
plt.plot(faxis, abs(np.fft.fft(filt_ecg)))
plt.xlim(-10,fs/2)
plt.xlabel('Frequency(Hz)')
plt.ylabel('')
plt.title('Filtered ECG data FFT')
plt.grid(which='major', color='black', linewidth=0.8)
plt.grid(which='minor', color='black', linewidth=0.5)
plt.minorticks_on()
plt.savefig('ecg_filtered_fft_%s.svg' %file, format = 'svg')

plt.show()