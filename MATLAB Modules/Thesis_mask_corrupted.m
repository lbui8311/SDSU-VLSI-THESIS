% Clear Workspace, commands and windows
clear; clc; close all;

%% 1. Reading BCI Data

% MEArec synthetic signal, background noise included
f  = 'C:\Users\Dr_ Luan\Documents\MATLAB\E218-BCI-Files\Sabes_Data\dataset_files\mearec_test_10s.h5';
X_in  = h5read(f, '/recordings');               
t = h5read(f, '/timestamps');                     % seconds, one value per sample
X_fs = double(h5read(f, '/info/recordings/fs'));   % 32 kHz
P  = h5read(f, '/channel_positions').';           % electrode x, y, z in um
st0 = h5read(f, '/spiketrains/0/times');          % neuron 0; folders are numbered 0, 1, 2, ...
duration = size(X_in,1) / X_fs;

% Update data retention to 60 seconds (matching numSamples = fs * 60)
fs = 10e3;
numSamples = fs * 10; % 10 instead of 60 for quick test purpose
numChannels = 10;
filtered_data = zeros(numChannels, numSamples);
X_v = X_in(1:numChannels,1:numSamples) * 1e-6; % Convert to volt
time_signal = (0:numSamples-1) / fs * 1e3;

for k = 1 : numChannels
    % Perform filtering step on all channels
    filtered_data(k, :) = FilterSBP(X_v(k, :), 10e3, 300, 3000, 'filtfilt', 'order', 2);
end

[s, f, t] = stft(X_in(1,:), fs, 'FrequencyRange', 'onesided');
[s_filt, f_filt, t_filt] = stft(filtered_data(1,:), fs, 'FrequencyRange', 'onesided');

% 2. Convert the complex STFT coefficients to dB magnitude
s_dB = 20*log10(abs(s));
s_dB_filt = 20*log10(abs(s_filt));

% =========================================================================
% FIND GLOBAL COLOR LIMITS (Forced to double to prevent type errors)
% =========================================================================
s_dB(isinf(s_dB)) = NaN;
s_dB_filt(isinf(s_dB_filt)) = NaN;

min_val = double(min([min(s_dB(:)), min(s_dB_filt(:))]));
max_val = double(max([max(s_dB(:)), max(s_dB_filt(:))]));
color_limits = [min_val, max_val];

%% 2. Masking the channels

% Channel Masking setup parameters
tau_sat = 1000 * 10^(-6); % Maximum Amplitude point (V)
var_min = 1e-11; % Lower variance threshold (V^2)
var_max = 1e-5; % Upper variance threshold (V^2)
snr_dB = 0; % Noise level: 0 dB = noise as strong as the signal
%snr_list = [0 -10 -20 -30 -40 -50]; % Noise level to test (dB)
snr_list = [0 -10 -20 -30 -40 -45 -50 -55 -60];
results = [];
noise_model = {'WGN', 'AR'};
win = round(0.1 * fs);  % mask window: 0.1 s = 1000 samples

X = filtered_data; % Clean signal

% [channel, first sample, last sample]
corrupt_events = [2 10001 12000;
                  4 30001 35000;
                  6 65001 69000;
                  8 90001 96000];

% Baseline: windowed mask on the clean signal (Amplitude Saturation + Variance Stability)
clean_flags = false(size(X));
% for s = 1:win:(size(X,2) - win + 1)
%     index = s : s + win - 1;
%     [~, ~, amplitudeMask, varianceMask, ~, clean_vf] = ...
%         Channel_Mask(X(:, index), tau_sat, var_min, var_max);
%     clean_flags(~(amplitudeMask & varianceMask), index) = true; % fails amplitude OR variance (1 = good)
% end


    [~, ~, amplitudeMask, varianceMask, ~, clean_vf] = ...
        Channel_Mask(X(:, :), tau_sat, var_min, var_max);
    %clean_flags(~(amplitudeMask & varianceMask), :) = true; % fails amplitude OR variance (1 = good)

% % warm-up: first sample where each channel becomes good; exclude it from scoring
% frameLen = size(X, 2) / size(clean_vf, 2);                  % = 100 samples = 10 ms
% clean_vf = repelem(clean_vf, 1, frameLen);                  % now 10 x 100000
% clean_mask = ~clean_vf;                                     % clean_vf: 1 = good  ->  clean_mask: true = bad
% 
% warm  = max(arrayfun(@(c) find(clean_vf(c, :), 1), 1:size(X, 1)));
% score = true(size(X));  score(:, 1:warm) = false;             % samples that count
% baseline_FMR = mean(clean_mask(score));

clean_flags = clean_vf;

baseline_FMR = mean(clean_flags(:)); % false-mask rate
fprintf('Baseline clean false-mask rate: %.4f\n', baseline_FMR);

for type = 1:length(noise_model)

    for snr_dB = snr_list

        % Noise injection
        rng(1); % Same noise every run, seed setup
        Xc = X; % Corrupted signal
        ground_mask = false(size(X)); % Ground truth flags: true = corrupted

        % Inject all events
        for k = 1: size(corrupt_events, 1)
            channel = corrupt_events(k, 1);
            start_sample = corrupt_events(k,2);
            end_sample = corrupt_events(k, 3);
            sample_length = end_sample - start_sample + 1;

            if type == 1
                % White Gaussian noise
                noiseSegment = randn(1, sample_length); 
            else
                % Auto-regressive noise
                noiseSegment = filter(1, [1 -0.95], randn(1, sample_length)); % y[k] = 0.95 * y[k-1] + e[k]
            end
            
            sigma = std(X(channel, :)) * 10^(-snr_dB / 20); % Auto-scale SNR to each channel's own signal
            noiseSegment = noiseSegment / std(noiseSegment) * sigma; % unit size

            Xc(channel, start_sample:end_sample) = Xc(channel, start_sample:end_sample) + noiseSegment; % add noise
            ground_mask(channel, start_sample:end_sample) = true;
        end

        % Run windowed mask on corrupted signal
        estimate_flags = false(size(Xc));
        % for s = 1:win:(size(Xc,2) - win + 1)
        %     index = s : s + win - 1;
        %     [~, ~, amplitudeMask, varianceMask, ~, vf] = ...
        %         Channel_Mask(Xc(:, index), tau_sat, var_min, var_max);
        %     estimate_flags(~(amplitudeMask & varianceMask), index) = true; % fails amplitude OR variance (1 = good)
        % end
            [~, ~, amplitudeMask, varianceMask, ~, vf] = ...
                Channel_Mask(Xc(:, :), tau_sat, var_min, var_max);
            estimate_flags(~(amplitudeMask & varianceMask), :) = true; % fails amplitude OR variance (1 = good)
   

        %estimate_mask = estimate_flags;
        estimate_mask = vf;

        % vf = repelem(vf, 1, frameLen);  % 10 x 100000
        % estimate_mask = ~vf; % vf: 1 = good  ->  estimate_mask: true = bad

        % Metrics 1-3
        true_negative = sum(~estimate_mask(:) & ~ground_mask(:)); % clean, kept
        false_negative = sum(~estimate_mask(:) & ground_mask(:)); % corrupted, missed
        false_positive = sum(estimate_mask(:) & ~ground_mask(:)); % clean, wrongly masked
        true_positive = sum(estimate_mask(:) & ground_mask(:));  % corrupted, masked 
 
        fault_recall = true_positive / (true_positive + false_negative); % of all corrupted samples, how many were caught
        mask_precision = true_positive / (true_positive + false_positive); % of all masked samples, how many were really corrupted
        clean_false_mask = false_positive / (false_positive + true_negative); % of all clean samples, how many were wrongly thrown away

        % Metrics 4-5
        onset = NaN(1, size(corrupt_events, 1));
        release = NaN(1, size(corrupt_events, 1));

        for k = 1:size(corrupt_events,1)
            eventChannel = corrupt_events(k, 1);
            eventStart = corrupt_events(k, 2);
            eventEnd = corrupt_events(k, 3);
            
            m = estimate_mask(eventChannel, :);
            i = find(m(eventStart:eventEnd), 1);

            if ~isempty(i)
                onset(k) = (i - 1) / fs;
                j = find(~m(eventEnd+1:end), 1); % First unmasked sample after it
                if ~isempty(j)
                    release(k) = (j - 1) / fs;
                else
                    release(k) = Inf;              % never released before the end of the file
                end
            end
        end

        onsetMs = NaN;  releaseMs = NaN;                  % NaN = nothing to time

        if any(~isnan(onset)),  onsetMs   = 1000 * median(onset(~isnan(onset)));       end
        if any(isfinite(release)), releaseMs = 1000 * median(release(isfinite(release))); end
        missed = sum(isnan(onset));  notReleased = sum(isinf(release));
        
        fprintf('%-3s %4d dB: recall %.3f | precision %.3f | clean false-mask %.4f | onset %6.1f ms | release %6.1f ms | missed %d | not released %d\n', ...
            noise_model{type}, snr_dB, fault_recall, mask_precision, clean_false_mask, onsetMs, releaseMs, missed, notReleased);

        results = [results; type, snr_dB, fault_recall, mask_precision, clean_false_mask, onsetMs, releaseMs, missed, notReleased];
 
    end
                                                                                                                                                                                                                                   
end

% grounded_FMR = mean(ground_mask(:)); % false-mask rate
% fprintf('Grounded false-mask rate: %.4f\n', grounded_FMR);
% 
% estimate_mask_FMR = mean(estimate_mask(:)); % false-mask rate
% fprintf('Corrupted false-mask rate: %.4f\n', estimate_mask_FMR);
% 
% [~, ~, ~, vr3, ~, ~] = Channel_Mask(X, tau_sat, var_min, 5e-11);   % 1e-11 < 5e-11 < clean variance
% vr3'     % all 0 -> upper limit acts on raw variance;  all 1 -> it doesn't
% 
% T = X;  T(7, :) = 0.03 * randn(1, size(X, 2));               % 30 mV noise, variance 9e-4
% [~, ~, am4, vr4, ~, ~] = Channel_Mask(T, tau_sat, var_min, var_max);
% [am4 vr4]     % row 7 should be 0 in vr4 if the upper variance check works

% V_clean = EMA_Running_Variance(filtered_data, 0.1);  
% fprintf('median clean variance: %.3g   (expect ~1.4e-10)\n', median(V_clean(:)));
% fprintf('below var_min: %.3f | above var_max: %.3f\n', ...
%     mean(V_clean(:) < var_min), mean(V_clean(:) > var_max));
% fprintf('fraction masked in clean_mask: %.3f\n', mean(clean_flags(:)));

% Save results table
T = array2table(results, 'VariableNames', {'model', 'snr_dB', 'fault_recall', 'mask_precision', ...
    'clean_false_mask', 'onset_ms', 'release_ms', 'missed', 'not_released'});
T.model = categorical(noise_model(T.model)');
writetable(T, 'mask_benchmark_results.csv');

% Plot recall vs SNR
figure; hold on
for type = 1:length(noise_model)
    sel = results(:, 1) == type;
    plot(results(sel, 2), results(sel, 3), '-o', 'DisplayName', noise_model{type});
end
set(gca, 'XDir', 'reverse'); grid on; legend
xlabel('SNR (dB)  (right = stronger noise)'); ylabel('Fault recall'); title('Detection vs noise level')

% Direct test to see if the mask detects flat channel and/or noisy channel
T = X;
T(3, :) = 0;                              % dead channel -> must fail the variance check (var < var_min)
T(7, :) = 5e-3 * randn(1, size(X, 2));    % 5 mV noise   -> must fail amplitude and variance
[~, m_all, am, vr, w] = Channel_Mask(T, tau_sat, var_min, var_max);

% % 1. Where is the mask on? (white = masked)
% figure; imagesc(estimate_mask); colormap(gray); xlabel('Sample'); ylabel('Channel')
% 
% % 2. How long is the start-up in this run vs. the clean baseline? (seconds)
% arrayfun(@(c) find(vf(c, :), 1), 1:size(X, 1)) / fs          % corrupted run
% arrayfun(@(c) find(clean_vf(c, :), 1), 1:size(X, 1)) / fs    % clean: ~0.6-0.75 s

%---------------------------------------------
% Single channel noise types comparison
channel = 1;
display_seconds = 0.05;                  % zoomed interval shown in the plot
num_samples = min(round(display_seconds * fs), size(X_in, 2));
clean_signal = X_in(channel, 1:num_samples) * 1e-6; % uV to V
time_ms = (0:num_samples-1) / fs * 1e3;

snr_dB = -10;                            % lower = stronger noise
rng(7);                                  % deterministic example

% Both noise types receive the same RMS level for a fair comparison.
noise_sigma = std(clean_signal) * 10^(-snr_dB / 20);
wgn = randn(size(clean_signal));
wgn = wgn / std(wgn) * noise_sigma;

ar_noise = filter(1, [1 -0.95], randn(size(clean_signal))); % AR(1)
ar_noise = ar_noise / std(ar_noise) * noise_sigma;

signal_wgn = clean_signal + wgn;
signal_ar = clean_signal + ar_noise;

all_signals_uV = [clean_signal, signal_wgn, signal_ar] * 1e6;
y_limits = [min(all_signals_uV), max(all_signals_uV)];
padding = 0.05 * max(diff(y_limits), eps);
y_limits = y_limits + [-padding, padding];

figure('Name', 'MEArec noise injection comparison', 'Color', 'w');
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(time_ms, clean_signal * 1e6, 'k', 'LineWidth', 1);
title('Clean MEArec signal');

nexttile;
plot(time_ms, signal_wgn * 1e6, 'b', 'LineWidth', 1);
title(sprintf('Gaussian noise (SNR = %d dB)', snr_dB));

nexttile;
plot(time_ms, signal_ar * 1e6, 'r', 'LineWidth', 1);
%title(sprintf('AR(1) noise, a = 0.95 (SNR = %d dB)', snr_dB));

for ax = findall(gcf, 'Type', 'axes').'
    ylim(ax, y_limits);
    xlabel(ax, 'Time (ms)');
    ylabel(ax, 'Voltage (uV)');
    grid(ax, 'on');
end

fprintf('Injected WGN and AR(1) noise at %d dB SNR (sigma = %.2f uV).\n', ...
    snr_dB, noise_sigma * 1e6);

%% Original vs Filtered (Voltage Amplitude & SNR baseline)

% First Subplot (Raw Data)
figure();
subplot(2,2,1);
plot(time_signal, X(1,:), '-k', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Voltage (V)');
title('Original MEArec Signal');
grid on
box on
set(gca, 'FontSize', 18, 'FontWeight', 'bold', 'FontAngle', 'italic');

subplot(2,2,2);
imagesc(t, f, s_dB);
axis xy; 
caxis(color_limits); 
xlabel('Time (s)');
ylabel('Frequency (Hz)');
set(gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');

% --- Single Colorbar Configuration ---
cb = colorbar; 
cb.Label.String = 'Power Magnitude (dB)'; 
cb.Label.FontSize = 18; 
cb.Label.FontWeight = 'bold';

% Second Subplot (Filtered Data)
subplot(2,2,3);
plot(time_signal, filtered_data(1,:), '-k', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Voltage (V)');
title('Filtered MEArec Signal');
grid on
box on
set(gca, 'FontSize', 18, 'FontWeight', 'bold', 'FontAngle', 'italic');

subplot(2,2,4);
imagesc(t_filt, f_filt, s_dB_filt);
axis xy; 
caxis(color_limits); 
xlabel('Time (s)');
ylabel('Frequency (Hz)');
set(gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');

% --- Single Colorbar Configuration ---
cb = colorbar; 
cb.Label.String = 'Power Magnitude (dB)'; 
cb.Label.FontSize = 18; 
cb.Label.FontWeight = 'bold';

% Original signal vs 10k SBP-filtered signal
% figure
% tiledlayout(1,1);
% nexttile;
% plot(time_signal, X(6,:), 'm', 'LineWidth', 2);
% hold on;
% plot(time_signal, filtered_data(6,:), 'c', 'LineWidth', 2);
% hold on;
% plot(time_signal, filtered_data_og(6,:), 'k', 'LineWidth', 2);
% hold off;
% xlabel('Time (s)'); ylabel('Voltage (V)');
% title(sprintf('Channel %d out of %d', 6, numChannels));
% legend('Original', 'SBP', 'SBP-OG');
% set(gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');
% grid on;