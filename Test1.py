"""
Created on Thu Oct 20 15:23:50 2022
@author: robert jack
"""

import numpy as np
import matplotlib.pyplot as plt


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

    def dofilter(self,ecgsig_val):
        self.buffer = np.roll(self.buffer,1)
        self.buffer[0] = ecgsig_val
        result = np.inner(self.buffer, self.coefficients)
        return result
            
    def doFilterAdaptive(self,signal,noise,learningRate):  #signal is the ecg data(the noisy signal)
        remover = self.dofilter(noise)
        output_e = signal - remover
        delta_h = output_e * learningRate * self.buffer
        self.coefficients = np.add(self.coefficients,delta_h)
        return output_e


data = np.loadtxt( 'robert_ecg.dat' )

print (data)
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
mu = 0.0005 #values larger than 0.001 made the system unstable
ntaps = 2000   #creating number of taps for lms filter
impulse_response, X = CreateImpulseResponse(1000, 45, 55, 2000)    
impulse_response_lms = np.zeros(ntaps)

filt_ecg = np.zeros(len(ecg_dat))
lms_filt_ecg = np.zeros(len(ecg_dat))

using_FIR = FIRfilter(impulse_response)  #instatiating the class

using_FIR_for_lms = FIRfilter(impulse_response_lms) 


inp_noise = np.zeros(len(ecg_dat))
for n in range(len(ecg_dat)):
    filt_ecg[n] = using_FIR.dofilter(ecg_dat[n]) #using filter to remove 50Hz, DC amd high frequencies
    
    inp_noise[n] = np.sin(2 * np.pi * n * fnoise/fs)
    lms_filt_ecg[n] = using_FIR_for_lms.doFilterAdaptive(ecg_dat[n], inp_noise[n], mu)
    inp_noise[n] = np.sin(2 * np.pi * n * fnoise2/fs)
    lms_filt_ecg[n] = using_FIR_for_lms.doFilterAdaptive(lms_filt_ecg[n], inp_noise[n], mu)
#cascading the lms filter to remove 50Hz and 210 Hz noise found from frequency spectrum 

#Creating the match filter
template = filt_ecg[1350:1500]
template_coeff = template[::-1]   

matchedFiltering = FIRfilter(template_coeff)
matchedFiltered = np.zeros(len(filt_ecg))
for i in range(len(filt_ecg)):
    matchedFiltered[i] = matchedFiltering.dofilter(filt_ecg[i])

matchedFiltered = matchedFiltered * matchedFiltered

peak = 0
last_peak = 0
beat_count = 0
inversebeatTime_array = np.zeros(len(filt_ecg))
for i in range(len(matchedFiltered)):
    if (matchedFiltered[i] >= 5) and (i - last_peak > (100*Ts)):
        beat_count += 1
        peak = i*Ts
        beatTime = peak - last_peak
        inverseBeatTime = 1/beatTime
        last_peak = peak
        inversebeatTime_array[i] = inverseBeatTime 



#################################PLOTTING####################################
    

fig1 = plt.figure(1)
plt.plot(faxis, abs(np.fft.fft(ecg_dat)))
plt.xlim(0,fs/2)
plt.xlabel('Frequency(Hz)')
plt.ylabel('')
plt.title('Initial ECG data FFT')
plt.grid(which='major', color='black', linewidth=0.8)
plt.grid(which='minor', color='black', linewidth=0.5)
plt.minorticks_on()
fig1.savefig('original_ecg_fft.pdf', format = 'pdf')

fig2 = plt.figure(2)
plt.plot(taxis, ecg_dat)
plt.xlabel('Time (s)')
plt.ylabel('Voltage (V)')
plt.title('Initial ECG data')
plt.grid(which='major', color='black', linewidth=0.8)
plt.grid(which='minor', color='black', linewidth=0.5)
plt.minorticks_on()
fig2.savefig('original_ecg.pdf', format = 'pdf')

fig3 = plt.figure(3)
plt.plot(impulse_response)
plt.xlabel('samples')
plt.ylabel('')
plt.title('Impulse Response')
plt.grid(which='major', color='black', linewidth=0.8)
plt.grid(which='minor', color='black', linewidth=0.5)
plt.minorticks_on()
fig3.savefig('impulse_response.pdf' , format = 'pdf')

fig4 = plt.figure(4)
plt.plot(X)
plt.xlabel('w (rad/s)')
plt.ylabel('')
plt.title('Bandstop')
plt.grid(which='major', color='black', linewidth=0.8)
plt.grid(which='minor', color='black', linewidth=0.5)
plt.minorticks_on()
fig4.savefig('Bandstop.pdf' , format = 'pdf')

fig5 = plt.figure(5)
plt.plot(taxis, filt_ecg)
plt.xlabel('Time(s)')
plt.ylabel('Voltage (V)')
plt.title('Filtered ECG data')
plt.grid(which='major', color='black', linewidth=0.8)
plt.grid(which='minor', color='black', linewidth=0.5)
plt.minorticks_on()
fig5.savefig('ecg_filtered.pdf' , format = 'pdf')

fig6 = plt.figure(6)
plt.plot(faxis, abs(np.fft.fft(filt_ecg)))
plt.xlim(0,fs/2)
plt.xlabel('Frequency(Hz)')
plt.ylabel('')
plt.title('Filtered ECG data FFT')
plt.grid(which='major', color='black', linewidth=0.8)
plt.grid(which='minor', color='black', linewidth=0.5)
plt.minorticks_on()
fig6.savefig('ecg_filtered_fft.pdf', format = 'pdf')


fig7 = plt.figure(7)
plt.plot(taxis, lms_filt_ecg) 
plt.xlabel('Time (s)')
plt.ylabel('Voltage (V)')
plt.title('lms Filtered ECG')
plt.grid(which='major', color='black', linewidth=0.8)
plt.grid(which='minor', color='black', linewidth=0.5)
plt.minorticks_on()
fig7.savefig('lms_filt_ecg.pdf' , format = 'pdf')

fig8 = plt.figure(8)
plt.plot(faxis, abs(np.fft.fft(lms_filt_ecg))) 
plt.xlabel('Frequency (Hz)')
plt.xlim(0,fs/2)
plt.ylabel('')
plt.title('lms Filtered ECG FFT')
plt.grid(which='major', color='black', linewidth=0.8)
plt.grid(which='minor', color='black', linewidth=0.5)
plt.minorticks_on()
fig8.savefig('lms_fft_ecg.pdf' , format = 'pdf')

fig9 = plt.figure(9)
plt.plot(taxis, matchedFiltered) 
plt.xlabel('Time (s)')
plt.ylabel('Voltage Squared (V**2)')
plt.title('Match Filtered ECG')
plt.grid(which='major', color='black', linewidth=0.8)
plt.grid(which='minor', color='black', linewidth=0.5)
plt.minorticks_on()
fig9.savefig('matchFiltered.pdf' , format = 'pdf')

fig10 = plt.figure(10)
plt.plot(faxis, abs(np.fft.fft(matchedFiltered))) 
plt.xlabel('Frequency (Hz)')
plt.xlim(0,fs/2)
plt.ylabel('')
plt.title('Match Filtered ECG FFT')
plt.grid(which='major', color='black', linewidth=0.8)
plt.grid(which='minor', color='black', linewidth=0.5)
plt.minorticks_on()
fig10.savefig('matchFiltered_fft_ecg.pdf' , format = 'pdf')

fig11 = plt.figure(11)
plt.plot(taxis, inversebeatTime_array) 
plt.xlabel('Time (s)')
plt.ylabel('')
plt.title('Inverse Beat Time ECG')
plt.grid(which='major', color='black', linewidth=0.8)
plt.grid(which='minor', color='black', linewidth=0.5)
plt.minorticks_on()
fig11.savefig('inverseBeatTime.pdf' , format = 'pdf')


plt.show()
