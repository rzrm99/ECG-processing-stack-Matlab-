# ECG-processing-stack

Full ECG processing stack: realistic signal synthesis, denoising, QRS detection, HRV metrics (time/frequency/nonlinear), a simple AF-like irregularity screen, and rich visualizations — all in `ecg_hrv_pipeline.m`.

---

# ECG + HRV Advanced Pipeline (MATLAB)

> **One-file demo** of a full ECG processing stack: realistic signal synthesis, denoising, QRS detection, HRV metrics (time/frequency/nonlinear), a simple AF-like irregularity screen, and rich visualizations — all in `ecg_hrv_pipeline.m`.

---

## 🧰 What’s inside

- **Signal synthesis**: Realistic ECG beats with RR variability (LF/HF) or AF-like jitter, plus common noise (baseline wander, powerline, white noise).
- **Denoising**: Baseline removal, 50/60 Hz notch, bandpass (QRS emphasis).
- **QRS detection**: Pan–Tompkins–style enhancement + adaptive thresholding and refractory handling; peak refinement on filtered ECG.
- **HRV analysis**:
  - Time domain: `meanNN`, `SDNN`, `RMSSD`, `pNN50`
  - Frequency domain (Welch PSD, LF/HF bands): `LFnu`, `HFnu`, band centers & powers
  - Nonlinear: Poincaré `SD1`, `SD2`
- **AF-like screen (heuristic)**: Flags when RR variability is high, autocorrelation low, and pNN50 elevated.
- **Visualization**: Noisy vs. clean ECG, filtered ECG with R-peaks, QRS-enhancement trace + threshold, tachogram, HRV PSD, Poincaré plot, LF/HF bars.

---

## ⚡ Quick start

1. Save the file as **`ecg_hrv_pipeline.m`**.
2. Open in **MATLAB** (R2019b+ recommended).
3. Press **Run** (F5).

The script detects whether the **Signal Processing Toolbox** is available.  
If yes, it uses `butter`, `filtfilt`, and `iirnotch`.  
If not, it falls back to robust moving-median/mean detrending and simple sinus fit for notch removal.

---

## ✅ Requirements

- **MATLAB** (tested with 360 Hz sampling; any modern version should work)
- **Optional**: Signal Processing Toolbox (for zero-phase IIR filtering)

No external data files needed — the ECG is synthesized.

---

## 🛠️ Usage

All parameters live at the top of the file:

```matlab
% -------------------- PARAMETERS --------------------
fs      = 360;      % sampling rate (Hz)
dur_s   = 180;      % duration (s)
hr_bpm  = 70;       % nominal heart rate
af_mode = false;    % true -> generate AF-like irregularity

% HRV (RR variability) settings
sdnn_ms = 60;       % target SDNN (ms)
lf_hz   = 0.1;      % LF modulation (~Mayer)
hf_hz   = 0.25;     % HF modulation (respiratory)

% Noise settings
add_baseline_wander = true;  baseline_amp_mV = 0.15; baseline_hz = 0.35;
add_powerline       = true;  powerline_freq  = 50;   powerline_amp_mV = 0.05;
```

---

## 📌 Examples

```text
Default (sinus rhythm)
    Just run ecg_hrv_pipeline.

AF-like rhythm
    Set af_mode = true.

North America mains
    Set powerline_freq = 60.

Shorter demo
    Set dur_s = 60.
```

---

## 📊 Outputs

**Printed in Command Window:**

```text
HRV metrics struct with fields:
meanNN_ms, sdnn_ms, rmssd_ms, pnn50, lfnu, hfnu,
lf_power, hf_power, total_power,
lf_center_hz, hf_center_hz,
SD1_ms, SD2_ms, cv_rr, ac_lag1, af_screen_positive
```

**Figures:**

1. Synthetic ECG with Noise (noisy + clean ground truth)  
2. Filtered ECG + R-peaks (stem markers)  
3. QRS Enhancement (integrated signal + threshold)  
4. Tachogram (RR in ms)  
5. HRV PSD (Welch) (0–0.5 Hz)  
6. Poincaré plot with SD1/SD2 annotation  
7. LF/HF normalized power bar chart  

---

## 🔎 Methods (high level)

```text
RR generation:
    Quasi-sinusoidal LF/HF modulations or AF-like random jitter;
    rescaled to target SDNN.

ECG morphology:
    Sum of Gaussian lobes for P/Q/R/S/T around R-times;
    slight beat-to-beat variability added.

Filtering:
    High-pass (~0.5 Hz) for baseline wander
    Notch at 50/60 Hz (Q≈30)
    Bandpass ~5–25 Hz for QRS emphasis

QRS detection:
    Differentiate → square → moving integration (~120 ms) →
    adaptive threshold + refractory (~200 ms) →
    local-max refinement on bandpassed ECG.

HRV frequency:
    Interpolate NN to 4 Hz; Welch PSD;
    integrate LF (0.04–0.15 Hz) and HF (0.15–0.40 Hz);
    compute normalized units and band centroids.

Poincaré:
    SD1/SD2 via rotated coordinates.

AF-like screen (heuristic):
    cv(RR), pNN50, and lag-1 autocorrelation thresholding
    (informal; not a medical diagnostic).
```

---

## 🔧 Key functions

```text
generate_rr_series      — builds RR series (LF/HF or AF-like)
synth_ecg_from_rr       — Gaussian-lobe beat template synthesis
hp_filter, notch_filter, bp_filter — denoising stages (with/without SPT)
detect_r_peaks, refine_r_peaks     — QRS enhancement + peak picking
hrv_freq, hrv_psd, centroid        — frequency-domain HRV utilities
poincare_SD1_SD2                   — nonlinear HRV metrics
autocorr_stat                      — lagged autocorrelation for RR series
ternary                            — tiny helper
```

---

## ⚠️ Notes & Caveats

```text
Educational/research demo only.
The AF-like screen is a simple heuristic and must NOT be used for clinical purposes.

Peak detection parameters (integration window, refractory time, bandpass edges)
are reasonable defaults for 360 Hz; adjust for other sampling rates or morphologies.

Frequency-domain metrics depend on duration; ≥3–5 minutes recommended for stable LF/HF estimates.

If you disable SPT, filtering still works but is cruder; expect slightly different peaks and HRV numbers.
```

---

## 🐛 Troubleshooting

```text
Too few NN intervals warning:
    Increase dur_s, reduce noise, or adjust filtering.

Missed/extra R-peaks:
    Tweak win_ms, refractory time in detect_r_peaks, or bandpass edges.

Powerline not fully removed:
    Ensure powerline_freq matches your region (50 vs 60 Hz).
    With SPT, notch quality is better.
```

---

## 📦 Repository layout (suggested)

```text
.
├── ecg_hrv_pipeline.m     # all code (top-level script + functions)
└── README.md              # this file
```

---

## 📜 License

Choose a license (e.g., MIT) and place it in `LICENSE`.

---

## ✨ Citation

If this helps your work, please cite the repository and acknowledge:  
**“ECG + HRV Advanced Pipeline (MATLAB) — synthetic ECG, denoising, QRS detection, and HRV metrics demo.”**

---

## 🚨 Medical Disclaimer


This software is provided strictly for educational and research purposes.
It is NOT intended for use in diagnosis, treatment, or any clinical decision-making.

The included AF-like irregularity screen is a simplified heuristic and should NOT be relied upon
for any medical conclusions.

Use at your own risk. The authors and contributors assume no responsibility
for any consequences arising from the use of this code.






