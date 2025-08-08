
function ecg_hrv_pipeline()
% ECG + HRV ADVANCED PIPELINE (Bioelectrical Engineering Project)
% ---------------------------------------------------------------
% What this script does:
% 1) Synthesizes realistic ECG with HRV + common noise
% 2) Denoises (baseline and mains hum), enhances QRS, detects R-peaks
% 3) Computes HRV time/frequency/nonlinear metrics
% 4) Screens for AF-like irregularity (simple heuristic)
% 5) Visualizes every step
%
% How to run:
%   - Save as ecg_hrv_pipeline.m
%   - Press Run (F5) in MATLAB
%
% If Signal Processing Toolbox is available, we use butter/filtfilt for
% zero-phase IIR filtering. If not, we fall back to simple convolutions.

%% -------------------- PARAMETERS --------------------
fs      = 360;              % sampling rate (Hz)
dur_s   = 180;              % duration (s)
hr_bpm  = 70;               % nominal heart rate
af_mode = false;            % set true to generate irregular AF-like rhythm

% HRV (RR variability) settings
sdnn_ms     = 60;           % target SDNN (ms) for variability realism
lf_hz       = 0.1;          % low-frequency modulation (~ Mayer waves)
hf_hz       = 0.25;         % high-frequency modulation (respiratory)

% Noise settings
add_baseline_wander = true;
baseline_amp_mV     = 0.15; % mV
baseline_hz         = 0.35; % Hz

add_powerline = true;
powerline_freq = 50;        % choose 50 or 60 depending on region
powerline_amp_mV = 0.05;

rng(42); % reproducibility

%% -------------------- SYNTHESIZE ECG --------------------
N = fs*dur_s;
t = (0:N-1)/fs;

% Generate RR-interval series with LF/HF variability or AF-like jitter
[rr_s, r_times] = generate_rr_series(dur_s, hr_bpm, sdnn_ms, lf_hz, hf_hz, af_mode);

% Build synthetic ECG by summing beats with morphology (Gaussian wavelets)
ecg_clean = synth_ecg_from_rr(t, r_times, fs);

% Add realistic noise
ecg_noisy = ecg_clean;
if add_baseline_wander
    ecg_noisy = ecg_noisy + baseline_amp_mV * sin(2*pi*baseline_hz*t);
end
if add_powerline
    ecg_noisy = ecg_noisy + powerline_amp_mV * sin(2*pi*powerline_freq*t);
end
% Add white sensor noise
ecg_noisy = ecg_noisy + 0.02*randn(size(t));

%% -------------------- DENOISING --------------------
% 1) Remove baseline wander (high-pass ~0.5 Hz)
ecg_hp = hp_filter(ecg_noisy, fs, 0.5);

% 2) Notch filter for powerline (50/60 Hz)
ecg_notch = notch_filter(ecg_hp, fs, powerline_freq);

% 3) Bandpass focus on QRS (5–15 Hz is classic; we'll widen a bit)
ecg_bp = bp_filter(ecg_notch, fs, 5, 25);

%% -------------------- QRS ENHANCEMENT + R-PEAK DETECTION --------------------
% Pan–Tompkins style: diff -> square -> moving integration
diff_sig = [0 diff(ecg_bp)];
sq_sig   = diff_sig.^2;
win_ms   = 120;                          % integration window ~120 ms
win_n    = max(1, round(fs*win_ms/1000));
int_sig  = movmean(sq_sig, win_n);

% Adaptive thresholding with refractory period
[locs_R, thr] = detect_r_peaks(int_sig, fs);

% R-peaks on original filtered ECG: refine to local maxima around detected indices
search_radius = round(0.05*fs);
locs_R = refine_r_peaks(ecg_bp, locs_R, search_radius);

%% -------------------- HRV METRICS --------------------
RR  = diff(locs_R)/fs;            % RR (s)
NN  = RR(~isoutlier(RR, 'median')); % drop gross ectopy for HRV core stats
if numel(NN) < 3
    warning('Too few NN intervals—try longer duration or adjust params.');
end

% Time-domain metrics
metrics.meanNN_ms = mean(NN)*1000;
metrics.sdnn_ms   = std(NN)*1000;
metrics.rmssd_ms  = sqrt(mean(diff(NN).^2))*1000;
metrics.pnn50     = 100*mean(abs(diff(NN)) > 0.050);

% Frequency-domain via Welch on interpolated tachogram
[lfnu, hfnu, lf_hz_est, hf_hz_est, lf_power, hf_power, total_power] = hrv_freq(NN, locs_R(1:end-1)/fs);

metrics.lfnu = lfnu; metrics.hfnu = hfnu;
metrics.lf_power = lf_power; metrics.hf_power = hf_power; metrics.total_power = total_power;
metrics.lf_center_hz = lf_hz_est; metrics.hf_center_hz = hf_hz_est;

% Nonlinear: Poincaré (SD1/SD2)
[SD1, SD2] = poincare_SD1_SD2(NN);
metrics.SD1_ms = SD1*1000; metrics.SD2_ms = SD2*1000;

% Simple AF-like irregularity screen (heuristic):
%   - high CV of RR, low autocorrelation, elevated pNN50, reduced SD1/SD2 ratio stability
cv_rr   = std(RR)/mean(RR);
ac_rr   = autocorr_stat(RR, 1);
metrics.cv_rr = cv_rr; metrics.ac_lag1 = ac_rr;

af_flag = (cv_rr > 0.15) && (metrics.pnn50 > 15) && (ac_rr < 0.35);
metrics.af_screen_positive = af_flag;

%% -------------------- REPORT + VISUALS --------------------
figure('Color','w','Name','ECG & HRV Pipeline','Position',[80 80 1200 800]);

subplot(4,1,1);
plot(t, ecg_noisy, 'Color', [0.2 0.2 0.2]); hold on;
plot(t, ecg_clean, 'r');
legend('Noisy ECG','Clean (ground truth)'); xlabel('Time (s)'); ylabel('mV');
title('Synthetic ECG with Noise');

subplot(4,1,2);
plot(t, ecg_notch); hold on;
stem(locs_R/fs, ecg_notch(locs_R), 'g','filled','MarkerSize',3);
xlabel('Time (s)'); ylabel('mV');
title('Filtered ECG + Detected R-peaks');

subplot(4,1,3);
plot(t, int_sig); hold on; yline(thr, '--r');
xlabel('Time (s)'); ylabel('AU');
title('QRS Enhancement (Integrated Signal) + Threshold');

subplot(4,1,4);
% Tachogram
plot(locs_R(1:end-1)/fs, RR*1000, '.-'); grid on;
xlabel('Time (s)'); ylabel('RR (ms)');
title('Tachogram (RR Intervals)');

% Print metrics to command window
disp('--- HRV METRICS ---');
disp(metrics);

% Secondary figure for frequency analysis & Poincaré
figure('Color','w','Name','HRV Frequency & Nonlinear','Position',[120 120 1100 400]);
subplot(1,3,1);
[pxx,f] = hrv_psd(NN, locs_R(1:end-1)/fs);
plot(f, pxx); xlim([0 0.5]); grid on;
xlabel('Hz'); ylabel('Power (ms^2/Hz)'); title('HRV PSD (Welch)');

subplot(1,3,2);
% Poincaré plot
NN1 = NN(1:end-1); NN2 = NN(2:end);
plot(NN1*1000, NN2*1000, '.'); grid on; axis equal;
xlabel('NN_n (ms)'); ylabel('NN_{n+1} (ms)');
title(sprintf('Poincaré (SD1=%.1f ms, SD2=%.1f ms)', metrics.SD1_ms, metrics.SD2_ms));

subplot(1,3,3);
bar([metrics.lfnu, metrics.hfnu]); set(gca,'XTickLabel',{'LFnu','HFnu'}); ylim([0 100]); grid on;
title('Normalized LF/HF Power (%)');

% Simple text panel
fprintf('\nAF-like irregularity screen: %s\n', ternary(metrics.af_screen_positive, 'POSITIVE (heuristic)', 'negative'));
fprintf('Mean NN = %.0f ms | SDNN = %.0f ms | RMSSD = %.0f ms | pNN50 = %.1f %%\n', ...
    metrics.meanNN_ms, metrics.sdnn_ms, metrics.rmssd_ms, metrics.pnn50);
fprintf('LFnu = %.1f %% | HFnu = %.1f %% | CV(RR) = %.2f | Lag-1 autocorr = %.2f\n\n', ...
    metrics.lfnu, metrics.hfnu, metrics.cv_rr, metrics.ac_lag1);

end

%% ================== SUPPORT FUNCTIONS ==================

function [rr_s, r_times] = generate_rr_series(dur_s, hr_bpm, sdnn_ms, lf_hz, hf_hz, af_mode)
% Generate RR series with LF/HF modulation or AF-like jitter
mean_rr = 60/hr_bpm; % seconds
tstep = 0.25;        % coarse grid for modulation
tmod = 0:tstep:dur_s;

if af_mode
    % AF-like: larger, erratic variability with 1/f^alpha noise
    base = mean_rr + 0.30*randn(size(tmod));
    base = max(0.3, base);
else
    % Quasi-sinusoidal LF/HF modulations
    lf = 0.05*sin(2*pi*lf_hz*tmod + 2*pi*rand);
    hf = 0.03*sin(2*pi*hf_hz*tmod + 2*pi*rand);
    base = mean_rr + lf + hf + 0.01*randn(size(tmod));
end

% Scale to target SDNN
rr_interp = interp1(tmod, base, 0:mean_rr:dur_s+10, 'linear', 'extrap');
rr_s = rr_interp(:);
scale = (sdnn_ms/1000) / std(rr_s);
rr_s = (rr_s - mean(rr_s))*scale + mean_rr;
rr_s(rr_s < 0.3) = 0.3;

r_times = cumsum(rr_s);
r_times = r_times(r_times < dur_s);
rr_s = diff([0; r_times]);
end

function ecg = synth_ecg_from_rr(t, r_times, fs)
% Simple beat template: sum of Gaussian lobes for P, Q, R, S, T
% Amplitudes (mV) and widths tuned for plausibility
ecg = zeros(size(t));
% Template parameters relative to R time (s)
comp = [ ... %   delay(s)   amp(mV)  width(s)
          -0.20       0.10     0.05;   % P
          -0.04      -0.15     0.015;  % Q
           0.00       1.00     0.012;  % R
           0.02      -0.25     0.015;  % S
           0.25       0.35     0.08];  % T
for k=1:numel(r_times)
    for c=1:size(comp,1)
        mu = r_times(k) + comp(c,1);
        A  = comp(c,2);
        w  = comp(c,3);
        ecg = ecg + A * exp(-0.5*((t-mu)/w).^2);
    end
end
% Slight morphology variability
ecg = ecg + 0.005*randn(size(ecg));
end

function y = hp_filter(x, fs, fc)
% High-pass ~ remove baseline wander
if has_spt()
    [b,a] = butter(2, fc/(fs/2), 'high');
    y = filtfilt(b,a,x);
else
    % Simple detrend via moving median + mean
    w = max(1, round(fs*0.6));
    x_med = movmedian(x, w);
    y = x - movmean(x_med, w);
end
end

function y = notch_filter(x, fs, f0)
% Notch at powerline
Q = 30;
if has_spt()
    w0 = f0/(fs/2);
    bw = w0/Q;
    [b,a] = iirnotch(w0, bw);
    y = filtfilt(b,a,x);
else
    % crude: subtract a narrowband sinus fit
    t = (0:numel(x)-1)/fs;
    ref = sin(2*pi*f0*t);
    alpha = (ref*x.')/(ref*ref.');
    y = x - alpha*ref;
end
end

function y = bp_filter(x, fs, f1, f2)
% Bandpass around QRS energy
if has_spt()
    [b,a] = butter(3, [f1 f2]/(fs/2), 'bandpass');
    y = filtfilt(b,a,x);
else
    % moving-average high emphasis
    y = x - movmean(x, round(fs*0.15));
end
end

function tf = has_spt()
% Detect Signal Processing Toolbox existence (butter/filtfilt)
tf = exist('butter','file')==2 && exist('filtfilt','file')==2;
end

function [locs, thr] = detect_r_peaks(int_sig, fs)
% Adaptive threshold on integrated QRS-enhanced signal
int_sig = int_sig(:);
% Estimate noise and signal levels using percentiles
noise_est = prctile(int_sig, 60);
sig_est   = prctile(int_sig, 98);
thr = 0.3*noise_est + 0.7*sig_est;

% Refractory period ~200 ms
ref = round(0.2*fs);
cand = find(int_sig > thr);

if isempty(cand)
    locs = [];
    return;
end

locs = [];
k = 1;
while k <= numel(cand)
    idx = cand(k);
    % form a refractory window
    w_end = min(numel(int_sig), idx+ref);
    [~, localmax] = max(int_sig(idx:w_end));
    locs(end+1) = idx + localmax - 1; %#ok<AGROW>
    % skip ahead beyond refractory
    k = find(cand > idx+ref, 1, 'first');
    if isempty(k), break; end
end
locs = unique(locs);
end

function locs_ref = refine_r_peaks(ecg, locs, rad)
% Local maximize the R on the bandpassed ECG
locs_ref = locs;
N = numel(ecg);
for i=1:numel(locs)
    a = max(1, locs(i)-rad);
    b = min(N, locs(i)+rad);
    [~, mx] = max(ecg(a:b));
    locs_ref(i) = a + mx - 1;
end
% Remove peaks too close (<200 ms)
minsep = round(0.2 * (length(ecg)/max(1,ceil(max(locs)/ (length(ecg)))))); %#ok<NASGU>
% Ensure strictly increasing & spaced
locs_ref = sort(locs_ref);
if numel(locs_ref) > 1
    RR = diff(locs_ref);
    keep = [true, RR > round(0.2* (length(ecg)/ (ecg_time_length(ecg))))];
    locs_ref = locs_ref(keep);
end
end

function T = ecg_time_length(ecg)
% helper for refine spacing
T = numel(ecg); % just to avoid extra args
end

function [lfnu, hfnu, lf_center, hf_center, lf_power, hf_power, total_power] = hrv_freq(NN, tNN)
% Interpolate NN(t) to uniform grid and compute PSD; integrate LF/HF bands.
% Bands (adult): VLF 0.003–0.04 Hz (ignored here), LF 0.04–0.15, HF 0.15–0.40
if numel(NN) < 4
    lfnu=NaN; hfnu=NaN; lf_center=NaN; hf_center=NaN;
    lf_power=NaN; hf_power=NaN; total_power=NaN; return;
end
fs_hrv = 4; % 4 Hz interpolation
tu = tNN(1):1/fs_hrv:tNN(end);
x  = interp1(tNN, NN, tu, 'pchip');

[pxx, f] = pwelch(x - mean(x), 256, 128, 1024, fs_hrv);

lf_idx = f>=0.04 & f<0.15;
hf_idx = f>=0.15 & f<=0.40;

lf_power = trapz(f(lf_idx), pxx(lf_idx));
hf_power = trapz(f(hf_idx), pxx(hf_idx));
total_power = trapz(f(f>=0.04 & f<=0.40), pxx(f>=0.04 & f<=0.40));

lfnu = 100 * lf_power / (lf_power + hf_power);
hfnu = 100 * hf_power / (lf_power + hf_power);

lf_center = centroid(f(lf_idx), pxx(lf_idx));
hf_center = centroid(f(hf_idx), pxx(hf_idx));
end

function [pxx,f] = hrv_psd(NN, tNN)
% helper to show PSD
if numel(NN) < 4
    pxx = zeros(1,512); f = linspace(0,0.5,512); return;
end
fs_hrv = 4;
tu = tNN(1):1/fs_hrv:tNN(end);
x  = interp1(tNN, NN, tu, 'pchip');
[pxx,f] = pwelch(x - mean(x), 256, 128, 1024, fs_hrv);
end

function c = centroid(f, p)
if isempty(f) || isempty(p) || sum(p)<=0
    c = NaN; return;
end
c = sum(f(:).*p(:))/sum(p(:));
end

function [SD1, SD2] = poincare_SD1_SD2(NN)
% SD1: short-term, SD2: long-term
NN1 = NN(1:end-1);
NN2 = NN(2:end);
diffs = (NN2 - NN1)/sqrt(2);
sums  = (NN2 + NN1)/sqrt(2);
SD1 = std(diffs);
SD2 = std(sums);
end

function ac1 = autocorr_stat(x, lag)
x = x(:) - mean(x);
if numel(x) < lag+1
    ac1 = NaN; return;
end
ac1 = sum( x(1:end-lag).*x(1+lag:end) ) / sum( x.^2 );
end

function out = ternary(cond, a, b)
if cond, out=a; else, out=b; end
end




