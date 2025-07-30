import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import welch
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from scipy.stats import entropy as scipy_entropy

"""
       Simulation of the two possibles patterns scenarios of the cell-cell interactions in pituitary, 
       reported in review: 
                    1) Signals showing only synchronicity with lower aperiodicity
                    2) Signals showing both asynchronicity and synchronicity (mixed interactions) between 
       Result plots where used for Figure 2 of  Review published in : 
"""

np.random.seed(42)
plt.style.use('fivethirtyeight')

# === Helper Functions ===

def generate_aperiodic_component(length, power_exp=1.5):
    white = np.random.randn(length)
    freqs = np.fft.rfftfreq(length)
    power = 1 / (freqs + 0.01)**power_exp
    fourier = np.fft.rfft(white) * power
    return np.fft.irfft(fourier, n=length)

def normalize(signal):
    return (signal - np.mean(signal)) / np.std(signal)

def generate_synchronous_signals(n_signals=3, length=2000, freq=10):
    t = np.linspace(0, 2, length)
    base_osc = np.sin(2 * np.pi * freq * t)
    signals = []
    for _ in range(n_signals):
        aperiodic = generate_aperiodic_component(length, power_exp=1.5)
        signal = normalize(base_osc + 0.4 * aperiodic)
        signals.append(signal)
    return np.array(signals)

def generate_mixed_signals(n_signals=3, length=2000, freq=10):
    t = np.linspace(0, 2, length)
    signals = []
    for i in range(n_signals):
        phase_shift = np.random.uniform(0, 2 * np.pi)
        freq_jitter = freq + np.random.randn() * 1.5
        amp_jitter = 0.4 + np.random.rand()
        osc = amp_jitter * np.sin(2 * np.pi * freq_jitter * t + phase_shift)
        noise_power = 0.5 if i % 2 == 0 else 1.5
        aperiodic = generate_aperiodic_component(length, power_exp=noise_power)
        signal = normalize(osc + 0.7 * aperiodic)
        signals.append(signal)
    return np.array(signals)

def phase_space_projection(signals):
    reshaped = signals.T
    pca = PCA(n_components=2)
    return pca.fit_transform(reshaped)

def compute_aperiodic_fit(signals, fs=1000):
    psds = []
    for signal in signals:
        freqs, psd = welch(signal, fs=fs, nperseg=512)
        psds.append(psd)
    psds = np.mean(psds, axis=0)

    log_freqs = np.log10(freqs[1:])
    log_psds = np.log10(psds[1:])

    model = LinearRegression()
    model.fit(log_freqs.reshape(-1, 1), log_psds)
    slope = model.coef_[0]

    return log_freqs, log_psds, slope

def compute_spectral_entropy(signals, fs=1000):
    """Average spectral entropy across all signals"""
    entropies = []
    for signal in signals:
        freqs, psd = welch(signal, fs=fs, nperseg=512)
        psd_norm = psd / np.sum(psd)
        ent = scipy_entropy(psd_norm)
        entropies.append(ent)
    return np.mean(entropies)

# === Simulation ===
length = 2000
fs = 1000

sync_signals = generate_synchronous_signals(length=length)
mixed_signals = generate_mixed_signals(length=length)

sync_corr = np.corrcoef(sync_signals)
mixed_corr = np.corrcoef(mixed_signals)

sync_proj = phase_space_projection(sync_signals)
mixed_proj = phase_space_projection(mixed_signals)

# Aperiodic slope
sync_freqs, sync_logpsd, sync_slope = compute_aperiodic_fit(sync_signals, fs=fs)
mixed_freqs, mixed_logpsd, mixed_slope = compute_aperiodic_fit(mixed_signals, fs=fs)

# Entropy
sync_entropy = compute_spectral_entropy(sync_signals, fs=fs)
mixed_entropy = compute_spectral_entropy(mixed_signals, fs=fs)

# === Plotting ===
fig, axs = plt.subplots(3, 2, figsize=(12, 14))

# Correlation Matrices
axs[0, 0].imshow(sync_corr, vmin=-1, vmax=1, cmap='coolwarm', alpha = 1.0)
axs[0, 0].set_title("Synchronous Correlation Matrix")
axs[0, 1].imshow(mixed_corr, vmin=-1, vmax=1, cmap='coolwarm')
axs[0, 1].set_title("Mixed Correlation Matrix")
fig.colorbar(axs[0, 0].images[0], ax=axs[0, :])

# Phase Space with Entropy annotation
axs[1, 0].plot(sync_proj[:, 0], sync_proj[:, 1], color='blue', lw=1.2)
axs[1, 0].set_title("Phase Space - Synchronous")
axs[1, 0].set_xlabel("PC1"); axs[1, 0].set_ylabel("PC2")
axs[1, 0].text(0.05, 0.95, f"Entropy: {sync_entropy:.3f}", transform=axs[1, 0].transAxes,
              fontsize=12, color='blue', ha='left', va='top', bbox=dict(boxstyle="round", fc="white"))

axs[1, 1].plot(mixed_proj[:, 0], mixed_proj[:, 1], color='red', lw=1.2)
axs[1, 1].set_title("Phase Space - Mixed Sync/Async")
axs[1, 1].set_xlabel("PC1"); axs[1, 1].set_ylabel("PC2")
axs[1, 1].text(0.05, 0.95, f"Entropy: {mixed_entropy:.3f}", transform=axs[1, 1].transAxes,
              fontsize=12, color='red', ha='left', va='top', bbox=dict(boxstyle="round", fc="white"))

# Aperiodic fits
axs[2, 0].plot(sync_freqs, sync_logpsd, label=f"Slope: {sync_slope:.2f}", color='blue')
axs[2, 0].set_title("Aperiodic Fit - Synchronous")
axs[2, 0].set_xlabel("log Frequency"); axs[2, 0].set_ylabel("log Power")
axs[2, 0].legend()

axs[2, 1].plot(mixed_freqs, mixed_logpsd, label=f"Slope: {mixed_slope:.2f}", color='red')
axs[2, 1].set_title("Aperiodic Fit - Mixed")
axs[2, 1].set_xlabel("log Frequency"); axs[2, 1].set_ylabel("log Power")
axs[2, 1].legend()

plt.tight_layout()
