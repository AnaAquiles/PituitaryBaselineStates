
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import spectrogram
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter
from scipy.signal import welch
from scipy.stats import norm


"""
      Simulations of aperiodic components and how different values are not showed in
      the periodic common representations, as 

       Result plots where used for Figure 3 and 4 of  Review published in : 
"""
plt.style.use('fivethirtyeight')

# Helper to create power-law noise (aperiodic 1/f^chi)
def generate_power_law_noise(exponent, size, fs):
    freqs = np.fft.rfftfreq(size, d=1/fs)
    freqs[0] = freqs[1]  # avoid division by zero
    amplitude = 1 / (freqs ** (exponent / 2))  # Power ~ 1/f^chi -> Amplitude ~ 1/f^(chi/2)
    phases = np.exp(1j * 2 * np.pi * np.random.rand(len(freqs)))
    spectrum = amplitude * phases
    signal = np.fft.irfft(spectrum, n=size)
    return signal

# Parameters
fs = 1000  # Hz
duration = 10  # seconds
n_samples = int(fs * duration)

# Aperiodic exponents (chi values)
chi_high = 1.2  # steeper, more low-frequency content
chi_low = 0.3   # flatter, more high-frequency content

# Generate signals
SignalHigh = []

for i in range(5):
    SignalHigh.append(generate_power_law_noise(chi_high, n_samples, fs))
SignalHigh = np.array(SignalHigh)  

SignalLow = []

for i in range(5):
    SignalLow.append(generate_power_law_noise(chi_low, n_samples, fs))
SignalLow = np.array(SignalLow) 


Spex = np.zeros((5001,len(SignalLow)))
Frex = np.zeros((5001,len(SignalLow)))

Cells = 5

plt.figure()
for d in range(Cells):
    # filter df for ticker and plot on specified axes
    Spex[:,d], Frex[:,d], _ = plt.magnitude_spectrum(SignalHigh[d,:], Fs = 2, color = 'slateblue', alpha = 0.6)
plt.xscale('log')
plt.yscale('log')
plt.title('High aperiodic value')

Spex_ = np.zeros((5001,len(SignalLow)))
Frex_ = np.zeros((5001,len(SignalLow)))

plt.figure()
for d in range(Cells):
    # filter df for ticker and plot on specified axes
    SpexL[:,d], FrexL[:,d], _ = plt.magnitude_spectrum(SignalLow[d,:], Fs = 2, color = 'slateblue', alpha = 0.2)
plt.xscale('log')
plt.yscale('log')
plt.title('Low aperiodic value')

## Selected whether activity you will measured High or Low aperiodic component 

Frecs = Frex   #FrexL for Low Aperiodic Value
Specs = Spex   #SpexL for Low Aperiodic Value 

fec = np.arange(0,1,0.01)
O_Freq = Frecs[:,0]
IndexF = np.array(np.where(O_Freq >= 0.01))
Index = IndexF.reshape(len(IndexF[0,:]))

Spex_ = Specs[Index]
Frex_ = Frecs[Index]

smoothPow = []

for i in range(Cells):
    smoothPow.append(gaussian_filter(Spex_[:,i], sigma = 2))

smoothPow = np.array(smoothPow)
powerMean = np.mean(smoothPow, axis = 0)
plt.figure()
for x in range(Cells):
    plt.plot(Frex_[:,x], smoothPow[x,:], c= 'gray', alpha =0.03, linewidth = 0.8)
    plt.plot(Frex_[:,0], powerMean, c ='r', linewidth = 1.5)
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Log power smoothed')

def combinedMod(f, alpha, k):
    return -(np.log(k + f**(alpha)))

# Fit the data
def fit_func(freqs, power):
    popt, pcov = curve_fit(combinedMod, freqs, power)
    alpha, k = popt
    return alpha

AperiodicExp = np.zeros(Cells)

for i in range(Cells):
    AperiodicExp[i] = np.abs(fit_func(Frex_[:,i], Spex_[:,i]))

"""
       Plot the transient signals simulated before helped with a butter bandpass filter
               going from 0.5 - 50 Hz , only for a best visualization 

"""

## Signal filtered 


fs = 1000
from scipy import signal
from scipy.signal import butter, lfilter

def butter_bandpass(lowcut, highcut, fs, order=5):
    return butter(order, [lowcut, highcut], fs=fs, btype='band')

def butter_bandpass_filter(data, lowcut, highcut, fs, order=3):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y

                        
lowcut = 0.5                                                        
highcut = 50

filtered = []
for i in range(0,5):
    filtered.append(butter_bandpass_filter(SignalLow[i], lowcut, highcut, fs, order = 3))
    
DataFiltBP = np.array(filtered)
 
# Time vector
t = np.linspace(0, duration, n_samples)

# Spectrograms
f_high, t_spec_high, Sxx_high = spectrogram(SignalHigh[0,:], fs=fs)
f_low, t_spec_low, Sxx_low = spectrogram(SignalLow[0,:], fs=fs)

# Plotting
fig, axs = plt.subplots(2, 2, figsize=(14, 8))

# Time series
axs[0, 0].plot(t, SignalHigh[0,:], color='red')
axs[0, 0].set_title('Signal with High Aperiodic Exponent (χ = 2.5)')
axs[0, 0].set_xlabel('Time (s)')
axs[0, 0].set_ylabel('Amplitude')

axs[1, 0].plot(t, SignalLow[0,:], color='blue')
axs[1, 0].set_title('Signal with Low Aperiodic Exponent (χ = 1.0)')
axs[1, 0].set_xlabel('Time (s)')
axs[1, 0].set_ylabel('Amplitude')

# Log-spectrograms
im1 = axs[0, 1].pcolormesh(t_spec_high, f_high, 10 * np.log10(Sxx_high + 1e-12), shading='gouraud')
axs[0, 1].set_title('Log-Spectrogram (χ = 2.5)')
axs[0, 1].set_xlabel('Time (s)')
axs[0, 1].set_ylabel('Frequency (Hz)')
fig.colorbar(im1, ax=axs[0, 1], label='Power (dB)')

im2 = axs[1, 1].pcolormesh(t_spec_low, f_low, 10 * np.log10(Sxx_low + 1e-12), shading='gouraud')
axs[1, 1].set_title('Log-Spectrogram (χ = 1.0)')
axs[1, 1].set_xlabel('Time (s)')
axs[1, 1].set_ylabel('Frequency (Hz)')
fig.colorbar(im2, ax=axs[1, 1], label='Power (dB)')

plt.tight_layout()
plt.show()

#### INDIVIDUAL PLOTS  (optional)
plt.figure()
plt.contourf(t_spec_high, f_high, 10 * np.log10(Sxx_high + 1e-12), shading='gouraud', cmap= 'jet')
plt.title('χ = 2.5')
plt.xlabel('Time (s)')
plt.ylabel('Frequency (Hz)')
plt.colorbar()

plt.figure()
plt.contourf(t_spec_low, f_low, 10 * np.log10(Sxx_low + 1e-12), shading='gouraud', cmap= 'jet')
plt.title('(χ = 1.0)')
plt.xlabel('Time (s)')
plt.ylabel('Frequency (Hz)')
plt.colorbar()

"""
          Toy example of the aperiodic components inspired in the figure described
      in  Donoghue, T., Haller, M., Peterson, E.J. et al.  DOI 10.1038/s41593-020-00744-x           

"""

# Parameters
fs = 1000  # Sampling frequency (Hz)
duration = 10  # seconds
n_samples = duration * fs
f_min, f_max = 1, 100  # Frequency range for analysis

# Generate a 1/f^β noise (aperiodic signal)
def generate_aperiodic_signal(beta=2.0):
    freqs = np.fft.rfftfreq(n_samples, 1/fs)
    amps = 1 / (freqs[1:] ** beta)
    phases = np.random.uniform(0, 2*np.pi, len(amps))
    spectrum = np.zeros(len(freqs), dtype=complex)
    spectrum[1:] = amps * np.exp(1j * phases)
    signal = np.fft.irfft(spectrum)
    return signal

# Generate signal
signal = generate_aperiodic_signal(beta=0.5)

# Compute PSD
freqs, psd = welch(signal, fs=fs, nperseg=2048)

# Only keep PSD in selected range
mask = (freqs >= f_min) & (freqs <= f_max)
freqs = freqs[mask]
psd = psd[mask]

# Fit log-log linear model: log10(PSD) = a * log10(f) + b
log_freqs = np.log10(freqs)
log_psd = np.log10(psd)
slope, offset = np.polyfit(log_freqs, log_psd, 1)

# Estimate -3dB bandwidth (relative to max power)
psd_max = np.max(psd)
threshold = psd_max / 2  # -3 dB
band_mask = psd >= threshold
bandwidth_freqs = freqs[band_mask]
bandwidth = bandwidth_freqs[-1] - bandwidth_freqs[0]

# Plot
plt.figure(figsize=(10, 6))
plt.loglog(freqs, psd, label='PSD')
plt.loglog(freqs, 10**(slope * log_freqs + offset), 'r--', label='Aperiodic Fit')

# Mark offset at 1 Hz
offset_freq = 1
offset_power = 10 ** (slope * np.log10(offset_freq) + offset)
plt.plot(offset_freq, offset_power, 'go', label=f'Offset (1 Hz): {offset_power:.2e}')

# Mark bandwidth range
plt.axvline(bandwidth_freqs[0], color='purple', linestyle=':', label='Bandwidth')
plt.axvline(bandwidth_freqs[-1], color='purple', linestyle=':')

plt.xlabel('Frequency (Hz)')
plt.ylabel('Power Spectral Density')
plt.title('Aperiodic Power Spectrum with Offset and Bandwidth')
plt.legend()
plt.grid(True, which="both", ls="--")
plt.tight_layout()
plt.show()

# Print results
print(f"Aperiodic offset (at 1 Hz): {offset_power:.3e}")
print(f"Slope (β): {-slope:.2f}")
print(f"Bandwidth (at -3 dB): {bandwidth:.2f} Hz")


"""


"""

plt.style.use('fivethirtyeight')

# Simulation parameters
fs = 1000
duration = 10
n_samples = fs * duration
freqs = np.fft.rfftfreq(n_samples, 1/fs)

# Create aperiodic (1/f^χ) component
def create_aperiodic_psd(frequencies, offset=1, exponent=1):
    return offset / (frequencies ** exponent)

# Add oscillatory peaks using Gaussians in frequency domain
def add_oscillatory_peaks(psd, frequencies, peaks):
    for cf, amp, bw in peaks:
        psd += amp * norm.pdf(frequencies, loc=cf, scale=bw/2.355)  # Convert FWHM to std
    return psd

# Generate frequency-domain PSD
aperiodic_exp = 1.5
offset = 1e-2
psd = create_aperiodic_psd(freqs[1:], offset=offset, exponent=aperiodic_exp)

# Define periodic peaks: (center freq, amplitude, FWHM)
peaks = [
    (10, 0.05, 4),   # Alpha peak
]
psd = add_oscillatory_peaks(psd, freqs[1:], peaks)

# Build complex spectrum with random phases
phases = np.exp(1j * 2 * np.pi * np.random.rand(len(psd)))
spectrum = np.zeros(len(freqs), dtype=complex)
spectrum[1:] = np.sqrt(psd) * phases
signal = np.fft.irfft(spectrum)

# Estimate PSD from time signal
f_welch, psd_welch = welch(signal, fs=fs, nperseg=2048)
log_f = np.log10(f_welch[1:])
log_psd = np.log10(psd_welch[1:])

# Fit linear model to log-log PSD (excluding peak region)
fit_mask = (f_welch > 2) & (f_welch < 40)
fit_mask &= ~((f_welch > 8) & (f_welch < 13))  # exclude alpha peak
slope, intercept = np.polyfit(np.log10(f_welch[fit_mask]), np.log10(psd_welch[fit_mask]), 1)

# Aperiodic model line
aperiodic_fit = 10 ** (intercept + slope * np.log10(f_welch))

# Plot
plt.figure(figsize=(10, 6))
plt.plot(f_welch, psd_welch, label="Empirical PSD", color='k')
plt.plot(f_welch, aperiodic_fit, 'r--', label=f"Aperiodic Fit\nSlope: {slope:.2f}, Offset: {10**intercept:.2e}")
plt.xscale("log")
plt.yscale("log")

# Highlight peak bandwidth (FWHM)
peak_cf, _, peak_bw = peaks[0]
plt.axvline(peak_cf - peak_bw / 2, color='blue', linestyle=':', label='Peak Bandwidth')
plt.axvline(peak_cf + peak_bw / 2, color='blue', linestyle=':')

# Offset marker at 1 Hz
offset_freq = 1
offset_power = 10 ** (intercept + slope * np.log10(offset_freq))
plt.plot(offset_freq, offset_power, 'go', label=f'Offset (1 Hz): {offset_power:.2e}')

plt.xlabel("Frequency (Hz)")
plt.ylabel("Power Spectral Density")
plt.title("Power Spectrum with Aperiodic Fit and Peak Bandwidth (Donoghue et al., 2020)")
plt.legend()
plt.grid(False)
plt.tight_layout()
plt.show()

"""

         The final function will simulate the signals again with different aperiodic values. 
         We will also add a correlation factor to represent diverse connectivity interactions,
                   which we will represent using common correlation matrices.

"""

# Helper to create power-law noise (aperiodic 1/f^chi)
def generate_power_law_noise(exponent, size, fs):
    freqs = np.fft.rfftfreq(size, d=1/fs)
    freqs[0] = freqs[1]  # avoid division by zero
    amplitude = 1 / (freqs ** (exponent / 2))  # Power ~ 1/f^chi -> Amplitude ~ 1/f^(chi/2)
    phases = np.exp(1j * 2 * np.pi * np.random.rand(len(freqs)))
    spectrum = amplitude * phases
    signal = np.fft.irfft(spectrum, n=size)
    return signal

# Parameters
fs = 1000  # Hz
duration = 10  # seconds
n_samples = int(fs * duration)

# Aperiodic exponents (chi values)
chi_high = 1.2  # steeper, more low-frequency content
chi_low = 0.3   # flatter, more high-frequency content

# Correlation factor
correlation_factor = 0.7  # Adjust this value to control the correlation between signals

# Generate signals
SignalHigh = []
for i in range(5):
    SignalHigh.append(generate_power_law_noise(chi_high, n_samples, fs))
SignalHigh = np.array(SignalHigh)

SignalLow = []
for i in range(5):
    SignalLow.append(generate_power_law_noise(chi_low, n_samples, fs))
SignalLow = np.array(SignalLow)

# Mix signals based on the correlation factor
for i in range(5):
    SignalLow[i] = correlation_factor * SignalHigh[i] + (1 - correlation_factor) * SignalLow[i]

# Compute correlation matrices
correlation_matrix_high = np.corrcoef(SignalHigh)
correlation_matrix_low = np.corrcoef(SignalLow)

# Plot correlation matrices
fig, axs = plt.subplots(1, 2, figsize=(16, 6))

# Plot for high signals
im_high = axs[0].imshow(correlation_matrix_high, cmap='jet', interpolation='nearest', alpha =0.5)
axs[0].set_title('Correlation Matrix of High Signals')
axs[0].set_xticks(np.arange(5))
axs[0].set_yticks(np.arange(5))
axs[0].set_xticklabels(np.arange(1, 6))
axs[0].set_yticklabels(np.arange(1, 6))
axs[0].set_xlabel('Signal Index')
axs[0].set_ylabel('Signal Index')
fig.colorbar(im_high, ax=axs[0], label='Correlation coefficient')

# Plot for low signals
im_low = axs[1].imshow(correlation_matrix_low, cmap='jet', interpolation='nearest', alpha =0.5)
axs[1].set_title('Correlation Matrix of Low Signals')
axs[1].set_xticks(np.arange(5))
axs[1].set_yticks(np.arange(5))
axs[1].set_xticklabels(np.arange(1, 6))
axs[1].set_yticklabels(np.arange(1, 6))
axs[1].set_xlabel('Signal Index')
axs[1].set_ylabel('Signal Index')
fig.colorbar(im_low, ax=axs[1], label='Correlation coefficient')

plt.tight_layout()
plt.show()


