% Clear Workspace, commands and windows
clear; clc; close all;

%% 1. Reading BCI Data

% MEArec synthetic signal, background noise included
f  = 'C:\Mac\Home\Documents\MATLAB\Thesis 2\E218-BCI-Files\Sabes_Data\dataset_files\mearec_test_10s.h5';
X_in  = h5read(f, '/recordings');               
t = h5read(f, '/timestamps');                     % seconds, one value per sample
X_fs = double(h5read(f, '/info/recordings/fs'));   % 32 kHz
P  = h5read(f, '/channel_positions').';           % electrode x, y, z in um
st0 = h5read(f, '/spiketrains/0/times');          % neuron 0; folders are numbered 0, 1, 2, ...

% % Sabes Data & Timestamps Paths
% data_file = 'C:\Users\Dr_ Luan\Documents\MATLAB\E218-BCI-Files\Sabes_Data\dataset_files\indy_20160624_03_fs10k.nwb';
% data_read_path = '/acquisition/timeseries/broadband/data';
% timestamps_read_path = '/acquisition/timeseries/broadband/timestamps';
% X_in = h5read(data_file, data_read_path);

duration = size(X_in, 2) / X_fs; 
fs = X_fs; % Keep 32 kHz sampling rate for accurate SBP filter
numSamples = min(size(X_in, 2), fs * 10); % 10 seconds of data
numChannels = 10;

filtered_data = zeros(numChannels, numSamples);
X_v = X_in(1:numChannels, 1:numSamples) * 1e-6; % Convert uV to V
time_signal = (0:numSamples-1) / fs; % time in seconds

for k = 1 : numChannels
    % Bandpass filter with actual sampling rate fs
    filtered_data(k, :) = FilterSBP(X_v(k, :), fs, 300, 3000, 'filtfilt', 'order', 2);
end

[s, f_stft, t_stft] = stft(X_v(1,:), fs, 'FrequencyRange', 'onesided');
[s_filt, f_filt, t_filt] = stft(filtered_data(1,:), fs, 'FrequencyRange', 'onesided');

s_dB = 20*log10(abs(s) + eps);
s_dB_filt = 20*log10(abs(s_filt) + eps);

s_dB(isinf(s_dB)) = NaN;
s_dB_filt(isinf(s_dB_filt)) = NaN;
min_val = double(min([min(s_dB(:)), min(s_dB_filt(:))]));
max_val = double(max([max(s_dB(:)), max(s_dB_filt(:))]));
color_limits = [min_val, max_val];

%% 2. Dynamic Threshold Estimation & Masking
X = filtered_data; % Clean baseline signal

% --- DYNAMIC THRESHOLD ESTIMATION ---
rms_per_channel = rms(X, 2);             % Baseline RMS per channel (V)
bg_rms = median(rms_per_channel);         % Median baseline RMS across channels (V)
bg_variance = bg_rms^2;                   % Baseline noise variance (V^2)

% Derive parameters dynamically from signal statistics
var_min = 0.15 * bg_variance;  % ~15% of baseline variance (dead channel cutoff)
var_max = 2.50 * bg_variance;  % ~2.5x baseline variance (noise burst cutoff)
tau_sat = 5.00 * bg_rms;       % ~5x baseline RMS voltage (saturation cutoff)

fprintf('--- Dynamically Estimated Thresholds ---\n');
fprintf('Median Baseline RMS : %6.2f uV\n', bg_rms * 1e6);
fprintf('Baseline Variance   : %.3e V^2\n', bg_variance);
fprintf('Estimated var_min   : %.3e V^2\n', var_min);
fprintf('Estimated var_max   : %.3e V^2\n', var_max);
fprintf('Estimated tau_sat   : %6.2f uV\n', tau_sat * 1e6);
fprintf('----------------------------------------\n\n');

snr_list = [20 15 10 6 3 0]; % 20 dB = subtle noise, 0 dB = strong noise
results = [];
noise_model = {'WGN', 'AR'};

% Corruption event definitions [channel, start_sample, end_sample]
corrupt_events = [2 10001 12000;
                  4 30001 35000;
                  6 65001 69000;
                  8 90001 96000];

% Baseline mask on clean data
[~, ~, ~, ~, ~, clean_vf] = Channel_Mask(X, tau_sat, var_min, var_max);
clean_flags = clean_vf; % Invert: true = good, false = corruped
baseline_FMR = mean(clean_flags(:));
fprintf('Baseline clean false-mask rate: %.4f\n', baseline_FMR);

for type = 1:length(noise_model)
    for snr_val = snr_list
        rng(1); 
        Xc = X; 
        ground_mask = false(size(X)); % Ground truth: true = corrupted

        % Inject noise events
        for k = 1:size(corrupt_events, 1)
            channel = corrupt_events(k, 1);
            start_sample = corrupt_events(k, 2);
            end_sample = corrupt_events(k, 3);
            sample_length = end_sample - start_sample + 1;

            if type == 1
                noiseSegment = randn(1, sample_length); 
            else
                noiseSegment = filter(1, [1 -0.95], randn(1, sample_length)); 
            end
            
            % Peak-to-peak scaled noise variance
            signal_p2p = peak2peak(X(channel, :));
            sigma = (signal_p2p / 6) * 10^(-snr_val / 20);
            noiseSegment = (noiseSegment / std(noiseSegment)) * sigma;

            Xc(channel, start_sample:end_sample) = Xc(channel, start_sample:end_sample) + noiseSegment;
            ground_mask(channel, start_sample:end_sample) = true;
        end

        % Run Channel Mask with estimated thresholds
        [~, ~, ~, ~, ~, vf] = Channel_Mask(Xc, tau_sat, var_min, var_max);
        estimate_mask = vf; % Invert validity flag (true = good, false =  bad/masked)

        % Metrics
        true_negative  = sum(~estimate_mask(:) & ~ground_mask(:)); % clean, kept
        false_negative = sum(~estimate_mask(:) &  ground_mask(:)); % corrupted, missed
        false_positive = sum( estimate_mask(:) & ~ground_mask(:)); % clean, wrongly masked
        true_positive  = sum( estimate_mask(:) &  ground_mask(:));  % corrupted, masked 
 
        fault_recall     = true_positive / max(1, (true_positive + false_negative)); % of all corrupted samples, how many were caught
        mask_precision   = true_positive / max(1, (true_positive + false_positive)); % of all masked samples, how many were really corrupted
        clean_false_mask = false_positive / max(1, (false_positive + true_negative)); % of all clean samples, how many were wrongly thrown away

        % Onset and release calculations
        onset = NaN(1, size(corrupt_events, 1));
        release = NaN(1, size(corrupt_events, 1));

        for k = 1:size(corrupt_events, 1)
            eventChannel = corrupt_events(k, 1);
            eventStart   = corrupt_events(k, 2);
            eventEnd     = corrupt_events(k, 3);
            
            m = estimate_mask(eventChannel, :);
            i = find(m(eventStart:eventEnd), 1);

            if ~isempty(i)
                onset(k) = (i - 1) / fs;
                j = find(~m(eventEnd+1:end), 1); 
                if ~isempty(j)
                    release(k) = (j - 1) / fs;
                else
                    release(k) = Inf; 
                end
            end
        end

        onsetMs = NaN; releaseMs = NaN;
        if any(~isnan(onset)),   onsetMs   = 1000 * median(onset(~isnan(onset))); end
        if any(isfinite(release)), releaseMs = 1000 * median(release(isfinite(release))); end
        missed = sum(isnan(onset)); 
        notReleased = sum(isinf(release));
        
        fprintf('%-3s %4d dB: recall %.3f | precision %.3f | clean false-mask %.4f | onset %6.1f ms | release %6.1f ms\n', ...
            noise_model{type}, snr_val, fault_recall, mask_precision, clean_false_mask, onsetMs, releaseMs);

        results = [results; type, snr_val, fault_recall, mask_precision, clean_false_mask, onsetMs, releaseMs, missed, notReleased];
    end
end

% Save results table
T = array2table(results, 'VariableNames', {'model', 'snr_dB', 'fault_recall', 'mask_precision', ...
    'clean_false_mask', 'onset_ms', 'release_ms', 'missed', 'not_released'});
T.model = categorical(noise_model(T.model)');
writetable(T, 'mask_benchmark_results.csv');

%% Original vs Filtered (Voltage Amplitude & SNR baseline)

% First Subplot (Raw Data)
figure();
subplot(2,2,1);
plot(time_signal, X_v(1,:), '-k', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Voltage (V)');
title('Original MEArec Signal');
grid on
box on
set(gca, 'FontSize', 18, 'FontWeight', 'bold', 'FontAngle', 'italic');

subplot(2,2,2);
imagesc(t_stft, f_stft, s_dB);
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
title('BPS-Filtered MEArec Signal');
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

%% Plot Recall vs SNR
figure('Color', 'w'); hold on;
for type = 1:length(noise_model)
    sel = results(:, 1) == type;
    plot(results(sel, 2), results(sel, 3), '-o', 'LineWidth', 1.5, 'DisplayName', noise_model{type});
end
set(gca, 'XDir', 'reverse'); % Left = High SNR (Clean), Right = Low SNR (Strong Noise)
grid on; legend('Location', 'best');
xlabel('SNR (dB) [Left = Clean (20dB), Right = Noisy (0dB)]', 'FontWeight', 'bold');
ylabel('Fault Recall', 'FontWeight', 'bold');
title('Detection Recall vs Noise Level', 'FontWeight', 'bold');

%% Multi-Channel Stacked Trace Plot
figure('Name', 'Multi-Channel Noise Injection Overlay', 'Color', 'w');
hold on;
time_sec = (0:numSamples - 1) / fs;
ch_spacing = max(range(Xc, 2)) * 1.2; 

for ch = 1:numChannels
    y_offset = (numChannels - ch) * ch_spacing;
    y_full = Xc(ch, :) + y_offset;
    
    plot(time_sec, y_full, 'k', 'LineWidth', 0.8);
    
    is_corrupt = ground_mask(ch, :);
    if any(is_corrupt)
        y_red = NaN(size(y_full));
        d = diff([0, is_corrupt, 0]);
        starts = find(d == 1);
        ends = find(d == -1) - 1;
        
        for e = 1:length(starts)
            idx_start = max(1, starts(e) - 1);
            idx_end = min(numSamples, ends(e) + 1);
            y_red(idx_start:idx_end) = y_full(idx_start:idx_end);
        end
        plot(time_sec, y_red, 'r', 'LineWidth', 1.2);
    end
end

xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Channel', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Multi-Channel Signal Overlay (%s Noise Model)', noise_model{type}), ...
    'FontSize', 13, 'FontWeight', 'bold');

yticks((0:numChannels - 1) * ch_spacing);
yticklabels(arrayfun(@(c) sprintf('Ch %d', c), numChannels:-1:1, 'UniformOutput', false));
xlim([time_sec(1), time_sec(end)]);
ylim([-ch_spacing * 0.5, numChannels * ch_spacing]);
grid on; box on;

h_clean = plot(NaN, NaN, 'k', 'LineWidth', 1);
h_noise = plot(NaN, NaN, 'r', 'LineWidth', 1.2);
legend([h_clean, h_noise], {'Clean Signal', 'Injected Noise'}, 'Location', 'northeast');



%% AR Noise modeling

figure ('Name', 'AR Noise')
tiledlayout(1,2);

% Tile 1: Original Signal
nexttile;
plot(time_signal, X_v(2,:), 'k', 'LineWidth', 2);
title('Unfiltered Original Signal')
xlabel('Time (s)'); ylabel('Voltage (V)');
set(gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');
grid on;

% Tile 2: Original Signal + AR Noise
nexttile;
plot(time_signal, Xc(2, :), 'k', 'LineWidth', 2);
title('Filtered Original Signal + AR Noise')
xlabel('Time (s)'); ylabel('Voltage (V)');
set(gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');
grid on;