% Clear Workspace, commands and windows
clear; clc; close all;

%% 2. Reading BCI Data
% Define Data & Timestamps Paths
data_file = 'C:\Users\Dr_ Luan\Documents\MATLAB\E218-BCI-Files\Sabes_Data\dataset_files\indy_20160624_03_fs10k.nwb';
data_read_path = '/acquisition/timeseries/broadband/data';
timestamps_read_path = '/acquisition/timeseries/broadband/timestamps';

% Read the Noisy_channels & Power_Noise_Channels
noisy_read = h5read(data_file, '/info/noisy_channels');
power_read = h5read(data_file, '/info/power_noise_channels');

% Update data retention to 60 seconds (matching numSamples = fs * 60)
fs = 10e3;
numSamples = fs * 10;

% Read the Data & Timestamps
voltage_signal = h5read(data_file, data_read_path);
voltage_signal = voltage_signal(1:2,1:numSamples);
time_signal = h5read(data_file, timestamps_read_path);
time_signal = time_signal(:,1:numSamples);

% Perform filtering step on all channels
numChannels = size(voltage_signal, 1);
filtered_data = zeros(numChannels, numSamples);

for k = 1 : numChannels
    filtered_data(k, :) = FilterSBP(voltage_signal(k, :), 10e3, 300, 3000, 'filtfilt', 'order', 2);
end

[s, f, t] = stft(voltage_signal(1,:), fs, 'FrequencyRange', 'onesided');
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

%% 3. FIR Filter: MATLAB vs. Verilog
% f_band = [8 60]; % Frequency bands range
% a = 1; % Denominator coefficients
% 
% % FIR 30th-order MATLAB generation
% coefficients_2 = fir1(30, f_band / (fs/2), 'bandpass');
% MAT_FIR = filter(coefficients_2, a, voltage_signal);
% 
%% 4. IIR Filter: MATLAB vs. Verilog
% % 6th-order Butterworth Filter (IIR) But Phase Distortion
% f_band = [8 60];
% [b, a] = butter(6, f_band / (fs/2), 'bandpass');
% 
% % 6th-order Butterworth Filter (IIR) using cascaded SOS (3)
% [sos, gain] = tf2sos(b, a); % SOS is 6 x 6
% 
% %[b01, b11, b21, a01, a11, a21]
% %[b02, b12, b22, a02, a12, a22]
% %[b03, b13, b23, a03, a13, a23]
% 
% % For causal, phase delay one-way filter: sosfilt()
% MAT_IIR = filtfilt(sos, gain, voltage_signal);


%% 3. Apply CAR AFTER filtering
data_mean = mean(filtered_data, 1);
data_car_all = filtered_data - data_mean;
data_car = data_car_all(1,:);

% % Artificial variance excursion
% data_car(10000:20000) = 0;
% data_car(35000:45000) = 0;
% data_car(50000:60000) = voltage_signal(50000:60000) * 50;
% 
% data_car(5000:7500) = 0;
% data_car(75000:85000) = 0;
% data_car(85000:90000) = 0;

%% 4. Amplitude Clipping

% Floating Point (64-bit Double)
capped_signal_float = zeros(numChannels, numSamples);

max_amp = 1000 * 10^(-6); % Maximum Amplitude point

for channel = 1:numChannels
    capped_signal_float = Amplitude_Clipping(data_car, max_amp);
end

% Fixed-Point (8/16/32-bit)
% Normalize to [-1,1] before Quantizing to avoid fractionals clip to zero
max_val = max(abs(data_car));
scale_factor = 1 / max_val;
scaled_signal = data_car * scale_factor;

max_amp_scaled = max_amp * scale_factor;

total_word_lengths = [8, 16, 32]; % Total bits
total_frac_lengths = [6, 14, 30]; % Fractional bits (2-bit for integers)
total_num_widths = length(total_word_lengths);

F = fimath(...
    'OverflowAction', 'Wrap', ...   % Allow overflow
    'RoundingMethod', 'Floor', ...  % Rounding down
    'SumMode', 'SpecifyPrecision', ...
    'ProductMode', 'SpecifyPrecision');

p = fipref;
p.FimathDisplay = 'none';
p.NumericTypeDisplay = 'none';

capped_signal_fixed = cell(total_num_widths, 1);
error_fixed = cell(total_num_widths, 1);
mse_fixed = zeros(total_num_widths, 1);
max_err_fixed = zeros(total_num_widths, 1);

fprintf('\nAmplitude Clipping errors:\n')

for k = 1:total_num_widths
    WL = total_word_lengths(k);
    FL = total_frac_lengths(k);

    F.SumWordLength = WL;
    F.ProductWordLength = WL;
    F.SumFractionLength = FL;
    F.ProductFractionLength = FL;

    % Quantize the signal
    signal_fixed = fi(scaled_signal, 1, WL, FL, 'fimath', F);

    % Run Amplitude Clipping
    signal_temp = Amplitude_Clipping(signal_fixed, max_amp_scaled);

    % Convert back to double
    capped_signal_fixed{k} = double(signal_temp) / scale_factor;

    assert(isequal(size(capped_signal_fixed{k}), size(capped_signal_float)), ...
    'Shape mismatch: clipped_fx{%d} is %s, clipped_f is %s', ...
    k, mat2str(size(capped_signal_fixed{k})), mat2str(size(capped_signal_float)));
    
    % Quantization error of fixed-point relative to floating-point
    error_fixed{k} = capped_signal_fixed{k} - capped_signal_float;
    mse_fixed(k) = mean(error_fixed{k}.^2);
    max_err_fixed(k) = max(abs(error_fixed{k}));

    fprintf('Q(%d,%d): MSE = %.4g V^2, Max Abs Error = %.4g V \n', ...
            WL, FL, mse_fixed(k), max_err_fixed(k));
end

% Plotting
colors = {'r', 'b', 'm'};
total_width_labels = arrayfun(@(w) sprintf('Fixed, %d-bit', w), total_word_lengths, 'UniformOutput', false);
% 
% figure ('Name', 'Amp Clip: Float vs Fixed')
% tiledlayout(2,2);
% 
% % Tile 1: Original Signal vs. Floating Point Clipped signals
% nexttile;
% plot(time_signal, data_car, 'b', 'LineWidth', 2);
% hold on;
% plot(time_signal, capped_signal_float, 'r', 'LineWidth', 2);
% hold off;
% title('Original vs. Floating Signals')
% xlabel('Time (s)'); ylabel('Voltage (V)');
% legend('Original', 'Clipped (Float, 64-bit)');
% set(gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');
% grid on;
% 
% % Tile 2: Floating vs Fixed Clipped signals overlaid
% nexttile;
% plot(time_signal, capped_signal_float, 'k', 'LineWidth', 2);
% hold on;
% for k = 1:total_num_widths
%     plot(time_signal, capped_signal_fixed{k}, colors{k}, 'LineWidth', 1);
% end
% hold off;
% title('Floating vs. Fixed Signals')
% xlabel('Time (s)'); ylabel('Voltage (V)');
% legend(['Float, 64-bit', total_width_labels]);
% set(gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');
% grid on;
% 
% % Tile 3: Quantization error for each word length
% nexttile;
% hold on;
% for k = 1:total_num_widths
%     plot(time_signal, error_fixed{k}, colors{k}, 'LineWidth', 1);
% end
% hold off;
% title('Quant Error (Max): Fixed vs Float')
% xlabel('Time (s)'); ylabel('Error (V)');
% legend(total_width_labels);
% set( gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');
% grid on;
% 
% % Tile 4: MSE vs. Total Word Lengths (Log-scale)
% nexttile;
% semilogy(total_word_lengths, mse_fixed, '-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
% title('Quant Error vs. Word Lengths')
% xlabel('Fixed-point Word Length (bits)'); ylabel('MSE (V^2)');
% set(gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');
% ylim([min(mse_fixed)*0.1, max(mse_fixed)*10]);
% grid on;

%% 5. SNR Proxy

signal_band = [300 3000]; % B_sig, 
noise_band = [300 3000]; % B_noise

% Floating Point (64-bit Double)
bandpower_float = zeros(numChannels, numSamples);
band_power_status_float = zeros(numChannels, 1);

for channel = 1:numChannels
    [band_power_status_float, bandpower_float] = SNR_Proxy(data_car, fs, signal_band, noise_band);
end

% Fixed-Point (8/16/32-bit)
bandpower_fixed = cell(total_num_widths, 1);
bandpower_status_fixed = cell(total_num_widths, 1);
error_fixed_SNR = cell(total_num_widths, 1);
mse_fixed_SNR = zeros(total_num_widths, 1);
max_err_fixed_SNR = zeros(total_num_widths, 1);

fprintf('\nSNR proxy errors:\n')

for k = 1:total_num_widths
    WL = total_word_lengths(k);
    FL = total_frac_lengths(k);

    F.SumWordLength = WL;
    F.ProductWordLength = WL;
    F.SumFractionLength = FL;
    F.ProductFractionLength = FL;

    % Quantize the signal
    signal_fixed_SNR = fi(scaled_signal, 1, WL, FL, 'fimath', F);

    % Run SNR Proxy
    [bandpower_status_fixed{k}, signal_temp_SNR] = SNR_Proxy(double(signal_fixed_SNR) / scale_factor, fs, signal_band, noise_band);

    n_frames = min(length(bandpower_float), length(signal_temp_SNR));

    % Convert back to double
    bandpower_fixed{k} = signal_temp_SNR(1:n_frames);

    assert(isequal(size(bandpower_fixed{k}), size(bandpower_float)), ...
        'Shape mismatch: clipped_fx{%d} is %s, clipped_f is %s', ...
        k, mat2str(size(bandpower_fixed{k})), mat2str(size(bandpower_float)));

    % Quantization error of fixed-point relative to floating-point
    error_fixed_SNR{k} = bandpower_fixed{k} - bandpower_float(1:n_frames);

    mse_fixed_SNR(k) = mean(error_fixed_SNR{k}.^2);
    max_err_fixed_SNR(k) = max(abs(error_fixed_SNR{k}));

    fprintf('Q(%d,%d): MSE = %.4g V^2, Max Abs Error = %.4g \n', ...
        WL, FL, mse_fixed_SNR(k), max_err_fixed_SNR(k));
end
% 
% figure ('Name', 'SNR Proxy: Float vs Fixed')
% tiledlayout(2,2);
% 
% % Tile 1: Floating Bandpower Ratio (SNR Proxy) per Frame
% nexttile;
% plot(1:length(bandpower_float), bandpower_float, 'r', 'LineWidth', 2);
% title('Floating SNR Proxy per Frame')
% xlabel('STFT Frame'); ylabel('SNR (dB)');
% yline(1.2, '--b', 'LineWidth',2); % tau_snr line
% legend('Float (64-bit)', '\tau_{snr} threshold');
% set(gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');
% grid on;
% 
% % Tile 2: Floating vs Fixed SNR Proxy
% nexttile;
% plot(1:length(bandpower_float), 10*log10(bandpower_float), 'k', 'LineWidth', 2);
% hold on;
% for k = 1:total_num_widths
%     plot(1:length(bandpower_fixed{k}), 10*log10(bandpower_fixed{k}), colors{k}, 'LineWidth', 1);
% end
% hold off;
% title('Floating vs. Fixed SNR Proxy')
% xlabel('STFT Frame'); ylabel('SNR (dB)');
% legend(['Float, 64-bit', total_width_labels]);
% set(gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');
% grid on;
% 
% % Tile 3: Quantization error for each word length
% nexttile;
% hold on;
% for k = 1:total_num_widths
%     plot(1:length(error_fixed_SNR{k}), error_fixed_SNR{k} , colors{k}, 'LineWidth', 1);
% end
% hold off;
% title('Quant Error (Max): Fixed vs Floating')
% xlabel('STFT Frame'); ylabel('Error (linear ratio)');
% legend(total_width_labels);
% set( gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');
% grid on;
% 
% % Tile 4: MSE vs. Total Word Lengths (Log-scale)
% nexttile;
% semilogy(total_word_lengths, mse_fixed_SNR, '-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
% title('SNR Quant Error vs. Word Lengths')
% xlabel('Fixed-point Word Length (bits)'); ylabel('MSE (V^2)');
% set(gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');
% ylim([min(mse_fixed_SNR)*0.1, max(mse_fixed_SNR)*10]);
% grid on;

%% 6. Variance Stability & Hysteretic Masking

% Function parameters
alpha = 0.1; % Smoothing factor for running mean and variance  
W = 10; % Sliding windows length in samples
max_fail = 5; % Consecutive fails before decay
max_pass = 10; % Consecutive pass before recovery
lambda = 0.9; % Weight decay factor if fail
delta = 0.05; % Weight recover factor if pass
% tau_min = 1e-13; % Sample Amplifier noise floor ~0.3 uV RMS (Dead/Shorted)
% tau_max = 1e-7; % Intra-cortical LFP rarely exceeds ~300 uV RMS (Massive Artifacts)

tau_min = 1e-11; 
tau_max = 1e-5; 

% Floating Point (64-bit Double)
VSHM_status_float = zeros(numChannels, 1);
VSHM_weights_float = zeros(numChannels, numSamples);

for channel = 1:numChannels
    [VSHM_status_float, VSHM_weights_float] = VS_HM(data_car, alpha, W, max_fail, max_pass, lambda, delta, tau_min, tau_max);
end

% Fixed-Point (8/16/32-bit)
VSHM_status_fixed = cell(total_num_widths, 1);
VSHM_weight_fixed = cell(total_num_widths, 1);
error_fixed_VSHM = cell(total_num_widths, 1);
mse_fixed_VSHM = zeros(total_num_widths, 1);
max_err_fixed_VSHM = zeros(total_num_widths, 1);

fprintf('\Variance Stability & Hysteretic Masking errors:\n')

for k = 1:total_num_widths
    WL = total_word_lengths(k);
    FL = total_frac_lengths(k);

    F.SumWordLength = WL;
    F.ProductWordLength = WL;
    F.SumFractionLength = FL;
    F.ProductFractionLength = FL;

    % Quantize the signal
    signal_fixed_VSHM = fi(scaled_signal, 1, WL, FL, 'fimath', F);

    % Run SNR Proxy
    [VSHM_status_fixed{k}, VSHM_weight_fixed{k}] = VS_HM(double(signal_fixed_VSHM) / scale_factor, alpha, W, max_fail, max_pass, lambda, delta, tau_min, tau_max);

    n_frames = min(length(VSHM_weights_float), length(VSHM_weight_fixed{k}));

    assert(isequal(size(VSHM_weight_fixed{k}), size(VSHM_weights_float)), ...
        'Shape mismatch: clipped_fx{%d} is %s, clipped_f is %s', ...
        k, mat2str(size(VSHM_weight_fixed{k})), mat2str(size(VSHM_weights_float)));

    % Quantization error of fixed-point relative to floating-point
    error_fixed_VSHM{k} = VSHM_weight_fixed{k}(1:n_frames) - VSHM_weights_float(1:n_frames);

    mse_fixed_VSHM(k) = mean(error_fixed_VSHM{k}.^2);
    max_err_fixed_VSHM(k) = max(abs(error_fixed_VSHM{k}));

    fprintf('Q(%d,%d): MSE = %.4g V^2, Max Abs Error = %.4g \n', ...
        WL, FL, mse_fixed_VSHM(k), max_err_fixed_VSHM(k));
end

figure ('Name', 'VSHM: Float vs Fixed')
tiledlayout(2,2);

% Tile 1: Floating channel weight trajectory
nexttile;
plot(1:length(VSHM_weights_float), VSHM_weights_float, 'r', 'LineWidth', 2);
title('Floating Channel Weight (w_i) per Sample')
xlabel('Sample'); ylabel('Weight w_i [0,1]');
ylim([-0.05 1.05]);
%legend('Float (64-bit)', '\tau_{snr} threshold');
set(gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');
grid on;

% Tile 2: Floating vs Fixed SNR Proxy
nexttile;
plot(1:length(VSHM_weights_float), VSHM_weights_float, 'k', 'LineWidth', 2);
hold on;
for k = 1:total_num_widths
    plot(1:length(VSHM_weight_fixed{k}),VSHM_weight_fixed{k}, colors{k}, 'LineWidth', 1);
end
hold off;
title('Channel Weight: Float vs Fixed')
xlabel('Sample'); ylabel('Weight w_i [0,1]');
legend(['Float, 64-bit', total_width_labels]);
ylim([-0.05 1.05]);
set(gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');
grid on;

% Tile 3: Quantization error for each word length
nexttile;
hold on;
for k = 1:total_num_widths
    plot(1:length(error_fixed_VSHM{k}), error_fixed_VSHM{k} , colors{k}, 'LineWidth', 1);
end
hold off;
title('Quant Error (Max): Fixed vs Float Weights')
xlabel('STFT Frame'); ylabel('Error (weight units)');
legend(total_width_labels);
set( gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');
grid on;

% Tile 4: MSE vs. Total Word Lengths (Log-scale)
nexttile;
%semilogy(total_word_lengths, mse_fixed_VSHM, '-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
bar(categorical(string(total_word_lengths) + "-bit"), mse_fixed_VSHM);
title('VSHM Quant Error vs. Word Lengths')
%xlabel('Fixed-point Word Length (bits)'); 
ylabel('MSE (V^2)');
set(gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');
if all(isfinite(mse_fixed_VSHM)) && max(mse_fixed_VSHM) > min(mse_fixed_VSHM)
    ylim([min(mse_fixed_VSHM)*0.1, max(mse_fixed_VSHM)*10]);
else
    warning('mse_fixed_VSHM has values: %s', mat2str(mse_fixed_VSHM));
end
grid on;

%% 7. Z-score normalization (Direct, Mitchell & PLA)

W = 10; % Windows size
alpha = 0.1; % Smoothing factor for running mean and variance
% epsilon = 1e-6; % Can't be 0 => inf
epsilon = 1e-14;

% Floating Point (64-bit Double)
z_direct_float = zeros(numChannels, numSamples);
z_mitchell_float = zeros(numChannels, numSamples);
z_pla_float = zeros(numChannels, numSamples);

for channel = 1:numChannels
    z_direct_float = Direct_Window_Z_Score(data_car, alpha, epsilon);
    z_mitchell_float = Mitchell_Approx_Z_Score(data_car, alpha, epsilon);
    [z_pla_float, seg_float] = PLA_LUT_Z_Score(data_car, alpha, epsilon);
end

[msg, warnID] = lastwarn;
disp(warnID);

% Fixed-Point (8/16/32-bit)
z_direct_fixed = cell(total_num_widths, 1);
z_mitchell_fixed = cell(total_num_widths, 1);
z_pla_fixed = cell(total_num_widths, 1);

% Need a settling period for EMA warmup instability 
warmup = 50;
error_fixed_direct_settled = cell(total_num_widths, 1);
error_fixed_mitchell_settled = cell(total_num_widths, 1);
error_fixed_pla_settled = cell(total_num_widths, 1);

error_fixed_direct = cell(total_num_widths, 1);
error_fixed_mitchell = cell(total_num_widths, 1);
error_fixed_pla = cell(total_num_widths, 1);

mse_fixed_direct = zeros(total_num_widths, 1);
mse_fixed_mitchell = zeros(total_num_widths, 1);
mse_fixed_pla = zeros(total_num_widths, 1);

max_err_fixed_direct = zeros(total_num_widths, 1);
max_err_fixed_mitchell = zeros(total_num_widths, 1);
max_err_fixed_pla = zeros(total_num_widths, 1);

% if k == 1
%     fprintf('\nZ-score Direct errors:\n')
% elseif k == 2
%     fprintf('\nZ-score Mitchell errors:\n')
% else
%     fprintf('\nZ-score PLA errors:\n')
% end

fprintf('\nZ-score Direct errors:\n')

seg_fixed = cell(total_num_widths, 1);
for k = 1:total_num_widths
    WL = total_word_lengths(k);
    FL = total_frac_lengths(k);

    F.SumWordLength = WL;
    F.ProductWordLength = WL;
    F.SumFractionLength = FL;
    F.ProductFractionLength = FL;

    % Quantize the signal
    signal_fixed_direct = fi(scaled_signal, 1, WL, FL, 'fimath', F);
    signal_fixed_mitchell = fi(scaled_signal, 1, WL, FL, 'fimath', F);
    signal_fixed_pla = fi(scaled_signal, 1, WL, FL, 'fimath', F);

    % Run Z-score
    z_direct_temp = Direct_Window_Z_Score(double(signal_fixed_direct) / scale_factor, alpha, epsilon);
    z_mitchell_temp = Mitchell_Approx_Z_Score(double(signal_fixed_mitchell) / scale_factor, alpha, epsilon);
    [z_pla_temp, seg_fixed{k}] = PLA_LUT_Z_Score(double(signal_fixed_pla) / scale_factor, alpha, epsilon);

    % Convert back to double
    % z_direct_fixed{k} = z_direct_temp / scale_factor;
    % z_mitchell_fixed{k} = z_mitchell_temp / scale_factor;
    % z_pla_fixed{k} = z_pla_temp / scale_factor;
    z_direct_fixed{k} = z_direct_temp;
    z_mitchell_fixed{k} = z_mitchell_temp;
    z_pla_fixed{k} = z_pla_temp;

    assert(isequal(size(z_direct_fixed{k}), size(z_direct_float)), ...
    'Shape mismatch: clipped_fx{%d} is %s, clipped_f is %s', ...
    k, mat2str(size(z_direct_fixed{k})), mat2str(size(z_direct_float)));

    assert(isequal(size(z_mitchell_fixed{k}), size(z_mitchell_float)), ...
        'Shape mismatch: clipped_fx{%d} is %s, clipped_f is %s', ...
        k, mat2str(size(z_mitchell_fixed{k})), mat2str(size(z_mitchell_float)));

    assert(isequal(size(z_pla_fixed{k}), size(z_pla_float)), ...
        'Shape mismatch: clipped_fx{%d} is %s, clipped_f is %s', ...
        k, mat2str(size(z_pla_fixed{k})), mat2str(size(z_pla_float)));
    
    % Quantization error of fixed-point relative to floating-point
    error_fixed_direct{k} = z_direct_fixed{k} - z_direct_float;
    error_fixed_direct_settled{k} = error_fixed_direct{k}(warmup+1:end);
    % mse_fixed_direct(k) = mean(error_fixed_direct{k}.^2);
    mse_fixed_direct(k) = mean(error_fixed_direct_settled{k}.^2);
    max_err_fixed_direct(k) = max(abs(error_fixed_direct_settled{k}));

    error_fixed_mitchell{k} = z_mitchell_fixed{k} - z_mitchell_float;
    error_fixed_mitchell_settled{k} = error_fixed_mitchell{k}(warmup+1:end);
    % mse_fixed_mitchell(k) = mean(error_fixed_mitchell{k}.^2);
    mse_fixed_mitchell(k) = mean(error_fixed_mitchell_settled{k}.^2);
    max_err_fixed_mitchell(k) = max(abs(error_fixed_mitchell_settled{k}));

    error_fixed_pla{k} = z_pla_fixed{k} - z_pla_float;
    error_fixed_pla_settled{k} = error_fixed_pla{k}(warmup+1:end);
    % mse_fixed_pla(k) = mean(error_fixed_pla{k}.^2);
    mse_fixed_pla(k) = mean(error_fixed_pla_settled{k}.^2);
    max_err_fixed_pla(k) = max(abs(error_fixed_pla_settled{k}));

    fprintf('Q(%d,%d) Direct: MSE = %.4g V^2, Max Abs Error = %.4g V \n', ...
            WL, FL, mse_fixed_direct(k), max_err_fixed_direct(k));

    fprintf('Q(%d,%d) Mitchell: MSE = %.4g V^2, Max Abs Error = %.4g V \n', ...
        WL, FL, mse_fixed_mitchell(k), max_err_fixed_mitchell(k));

    fprintf('Q(%d,%d) PLA: MSE = %.4g V^2, Max Abs Error = %.4g V \n', ...
        WL, FL, mse_fixed_pla(k), max_err_fixed_pla(k));

end

warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');

% Verifying  
segment_difference = (seg_float ~= seg_fixed{2});
fprintf('\n\nTotal samples crossed segment:');
disp(sum(segment_difference)); % How many samples crossed the segment boundary?
[~, spike_index] = findpeaks(abs(error_fixed_pla{2}), 'MinPeakHeight', 0.01); % Locate the visible spikes

fprintf('Total samples crossed peak 0.01:');
%disp(segment_difference(spike_index)); % Were those spikes among the sample crossed?
disp(length(segment_difference(spike_index)));

fprintf('Indices of those crossed segment:');
disp(find(segment_difference));           % all 12 sample indices where segment changed

fprintf('Error magnitudes of those crossed segment:');
disp(error_fixed_pla{2}(find(segment_difference)));   % actual error magnitude at all 12, not just the 10 spikes already found

% 
% figure ('Name', 'Z-score: Float vs Fixed (1)')
% tiledlayout(2,2);
% 
% % Tile 1: Original Signal vs. Floating Point Clipped signals
% nexttile;
% plot(time_signal, data_car, 'k', 'LineWidth', 2);
% hold on;
% plot(time_signal, z_direct_float, 'r', 'LineWidth', 2);
% hold on;
% plot(time_signal, z_mitchell_float, 'b', 'LineWidth', 2);
% hold on;
% plot(time_signal, z_pla_float, 'g', 'LineWidth', 2);
% hold off;
% title('Original vs. 64-bit Floating Signals')
% xlabel('Time (s)'); ylabel('Voltage (V)');
% legend('Original', 'Direct', 'Mitchell', 'PLA');
% set(gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');
% grid on;
% 
% % Tile 2: Direct Float vs Fixed signals overlaid
% nexttile;
% plot(time_signal, z_direct_float, 'k', 'LineWidth', 2);
% hold on;
% for k = 1:total_num_widths
%     plot(time_signal, z_direct_fixed{k}, colors{k}, 'LineWidth', 1);
% end
% hold off;
% title('Direct Float vs. Fixed Signals')
% xlabel('Time (s)'); ylabel('Voltage (V)');
% legend(['Float, 64-bit', total_width_labels]);
% set(gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');
% grid on;
% 
% % Tile 3: Quantization error for each word length
% nexttile;
% hold on;
% for k = 1:total_num_widths
%     plot(time_signal, error_fixed_direct{k}, colors{k}, 'LineWidth', 1);
% end
% hold off;
% title('Direct Quant Error (Max): Fixed vs Float')
% xlabel('Time (s)'); ylabel('Error (V)');
% legend(total_width_labels);
% set( gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');
% grid on;
% 
% % Tile 4: MSE vs. Total Word Lengths (Log-scale)
% nexttile;
% semilogy(total_word_lengths, mse_fixed_direct, '-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
% title('Direct Quant Error vs. Word Lengths')
% xlabel('Fixed-point Word Length (bits)'); ylabel('MSE (V^2)');
% set(gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');
% ylim([min(mse_fixed_direct)*0.1, max(mse_fixed_direct)*10]);
% grid on;
% 
% figure ('Name', 'Z-score: Float vs Fixed (2)')
% tiledlayout(2,2);
% 
% % Tile 1: Original Signal vs. Floating Point Clipped signals
% nexttile;
% plot(time_signal, data_car, 'k', 'LineWidth', 2);
% hold on;
% plot(time_signal, z_direct_float, 'r', 'LineWidth', 2);
% hold on;
% plot(time_signal, z_mitchell_float, 'b', 'LineWidth', 2);
% hold on;
% plot(time_signal, z_pla_float, 'g', 'LineWidth', 2);
% hold off;
% title('Original vs. 64-bit Floating Signals')
% xlabel('Time (s)'); ylabel('Voltage (V)');
% legend('Original', 'Direct', 'Mitchell', 'PLA');
% set(gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');
% grid on;
% 
% % Tile 2: Mitchell Float vs Fixed signals overlaid
% nexttile;
% plot(time_signal, z_mitchell_float, 'k', 'LineWidth', 2);
% hold on;
% for k = 1:total_num_widths
%     plot(time_signal, z_mitchell_fixed{k}, colors{k}, 'LineWidth', 1);
% end
% hold off;
% title('Mitchell Float vs. Fixed Signals')
% xlabel('Time (s)'); ylabel('Voltage (V)');
% legend(['Float, 64-bit', total_width_labels]);
% set(gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');
% grid on;
% 
% % Tile 3: Quantization error for each word length
% nexttile;
% hold on;
% for k = 1:total_num_widths
%     plot(time_signal, error_fixed_mitchell{k}, colors{k}, 'LineWidth', 1);
% end
% hold off;
% title('Mitchell Quant Error (Max): Fixed vs Float')
% xlabel('Time (s)'); ylabel('Error (V)');
% legend(total_width_labels);
% set( gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');
% grid on;
% 
% % Tile 4: MSE vs. Total Word Lengths (Log-scale)
% nexttile;
% semilogy(total_word_lengths, mse_fixed_mitchell, '-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
% title('Mitchell Quant Error vs. Word Lengths')
% xlabel('Fixed-point Word Length (bits)'); ylabel('MSE (V^2)');
% set(gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');
% ylim([min(mse_fixed_mitchell)*0.1, max(mse_fixed_mitchell)*10]);
% grid on;
% 
% figure ('Name', 'Z-score: Float vs Fixed (3)')
% tiledlayout(2,2);
% 
% % Tile 1: Original Signal vs. Floating Point Clipped signals
% nexttile;
% plot(time_signal, data_car, 'k', 'LineWidth', 2);
% hold on;
% plot(time_signal, z_direct_float, 'r', 'LineWidth', 2);
% hold on;
% plot(time_signal, z_mitchell_float, 'b', 'LineWidth', 2);
% hold on;
% plot(time_signal, z_pla_float, 'g', 'LineWidth', 2);
% hold off;
% title('Original vs. 64-bit Floating Signals')
% xlabel('Time (s)'); ylabel('Voltage (V)');
% legend('Original', 'Direct', 'Mitchell', 'PLA');
% set(gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');
% grid on;
% 
% 
% % Tile 2: PLA Floating vs Fixed Clipped signals overlaid
% nexttile;
% plot(time_signal, z_pla_float, 'k', 'LineWidth', 2);
% hold on;
% for k = 1:total_num_widths
%     plot(time_signal, z_pla_fixed{k}, colors{k}, 'LineWidth', 1);
% end
% hold off;
% title('PLA Float vs. Fixed Signals')
% xlabel('Time (s)'); ylabel('Voltage (V)');
% legend(['Float, 64-bit', total_width_labels]);
% set(gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');
% grid on;
% 
% 
% % Tile 3: Quantization error for each word length
% nexttile;
% hold on;
% for k = 1:total_num_widths
%     plot(time_signal, error_fixed_pla{k}, colors{k}, 'LineWidth', 1);
% end
% hold off;
% title('PLA Quant Error (Max): Fixed vs Float')
% xlabel('Time (s)'); ylabel('Error (V)');
% legend(total_width_labels);
% set( gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');
% grid on;
% 
% % Tile 4: MSE vs. Total Word Lengths (Log-scale)
% nexttile;
% semilogy(total_word_lengths, mse_fixed_pla, '-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
% title('PLA Quant Error vs. Word Lengths')
% xlabel('Fixed-point Word Length (bits)'); ylabel('MSE (V^2)');
% set(gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');
% ylim([min(mse_fixed_pla)*0.1, max(mse_fixed_pla)*10]);
% grid on;
% % 
% % disp(isequal(z_mitchell_fixed{1}, z_mitchell_fixed{3}));
% % 
% % disp(mean(abs(z_pla_float - z_direct_float)));   % PLA vs. true z-score, not vs. its own fixed-point
% % 
% % [~, idx] = max(abs(error_fixed_direct{1}));
% % disp(idx);                          % where is the single largest error?
% % disp(sum(abs(error_fixed_direct{1}) > 1e-4));   % how many samples exceed a "should be quantization-only" threshold?
% % 
% % disp(z_direct_float(1) - z_direct_fixed{1}(1));
% % disp(z_direct_float(1) - z_direct_fixed{3}(1));
% 
% format long
% disp([mse_fixed_mitchell(1), mse_fixed_mitchell(2), mse_fixed_mitchell(3)]);
% format short
% 
% disp(z_direct_float(1:5));   % if you can access this, or recompute it directly
% 
% disp(mse_fixed_direct);           % should now differ from your original flat 3.856e-05 result
% disp(mean(abs(z_pla_float - z_direct_float)));   % should shrink from 0.0041 if the segment_idx fix was applied
% 
% figure('Name', 'Direct: Float vs. Each Fixed-Point Width');
% tiledlayout(2,2);
% 
% nexttile;
% plot(time_signal, z_direct_float, 'k', 'LineWidth', 1);
% title('Floating-point (64-bit)');
% ylim([-4 4]); grid on;
% 
% for k = 1:total_num_widths
%     nexttile;
%     plot(time_signal, z_direct_float, 'k', 'LineWidth', 0.5); hold on;
%     plot(time_signal, z_direct_fixed{k}, colors{k}, 'LineWidth', 1); hold off;
%     title(sprintf('%d-bit vs. Float', total_word_lengths(k)));
%     ylim([-4 4]); grid on;
% end
% 
% warning('off', 'MATLAB:polyfit:PolyNotUnique');
% disp('Test below:');
% disp([mse_fixed_direct, mse_fixed_mitchell, mse_fixed_pla]);
% disp(mean(abs(z_pla_float - z_direct_float)));
% 
% disp(mean(abs(z_pla_float - z_direct_float)));  % overall, already have this
% figure; plot(abs(z_pla_float - z_direct_float)); title('PLA vs Direct disagreement over time');
% 
% figure('Name', 'Quantization Error Comparison Across Methods');
% semilogy(total_word_lengths, mse_fixed_direct, '-o', 'LineWidth', 2, 'Color', 'r'); hold on;
% semilogy(total_word_lengths, mse_fixed_mitchell, '-s', 'LineWidth', 2, 'Color', 'b');
% semilogy(total_word_lengths, mse_fixed_pla, '-^', 'LineWidth', 2, 'Color', 'g');
% hold off;
% xlabel('Fixed-point word length (bits)');
% ylabel('MSE (V^2, log scale)');
% title('Fixed-point quantization error by method');
% legend('Direct', 'Mitchell', 'PLA', 'Location', 'best');
% grid on;
% 
% accuracy_gap = [mean(abs(z_mitchell_float - z_direct_float)), ...
%     mean(abs(z_pla_float - z_direct_float))];
% 
% figure('Name', 'Floating-point approximation accuracy');
% bar(categorical({'Mitchell','PLA'}), accuracy_gap);
% ylabel('Mean |error| vs. direct-window ground truth (V)');
% title('Inherent approximation error, independent of quantization');
% grid on;
% 
% mse_matrix = [mse_fixed_direct'; mse_fixed_mitchell'; mse_fixed_pla'];  % 3 methods x 3 widths
% 
% figure('Name', 'MSE by method and word length');
% b = bar(categorical({'8-bit','16-bit','32-bit'}), mse_matrix', 'grouped');
% set(gca, 'YScale', 'log');
% ylabel('MSE (V^2, log scale)');
% legend('Direct','Mitchell','PLA', 'Location', 'best');
% title('Quantization error grouped by word length');
% grid on;
% 
% % Add value labels above each bar
% for i = 1:numel(b)
%     xtips = b(i).XEndPoints;
%     ytips = b(i).YEndPoints;
%     labels = arrayfun(@(v) sprintf('%.2g', v), b(i).YData, 'UniformOutput', false);
%     text(xtips, ytips, labels, 'HorizontalAlignment', 'center', ...
%     'VerticalAlignment', 'bottom', 'FontSize', 8, 'FontWeight', 'bold');
% end

%% Figure 1: Overview -- original signal and all three float z-scores
% Split into two panels since data_car (raw V) and the z-scores
% (dimensionless, ~0.8-3.2) are not on comparable scales.
figure('Name', 'Z-score: Overview');
tiledlayout(2,1);
 
nexttile;
plot(time_signal, data_car, 'k', 'LineWidth', 1);
title('Original signal (64-bit float)');
xlabel('Time (s)'); ylabel('Voltage (V)');
grid on;
 
nexttile;
plot(time_signal, z_direct_float, 'r', 'LineWidth', 1); hold on;
plot(time_signal, z_mitchell_float, 'b', 'LineWidth', 1);
plot(time_signal, z_pla_float, 'g', 'LineWidth', 1);
hold off;
title('Z-score, all methods (64-bit float)');
xlabel('Time (s)'); ylabel('Z-score (dimensionless)');
legend('Direct', 'Mitchell', 'PLA', 'Location', 'best');
grid on;
 
%% Figures 2-4: Per-method small multiples -- float vs. each word length
% One figure per method, 2x2: float alone, then float vs. 8/16/32-bit.
method_names   = {'Direct', 'Mitchell', 'PLA'};
z_float_all    = {z_direct_float, z_mitchell_float, z_pla_float};
z_fixed_all    = {z_direct_fixed, z_mitchell_fixed, z_pla_fixed};
 
for m = 1:3
    figure('Name', sprintf('%s: Float vs. Each Fixed-Point Width', method_names{m}));
    tiledlayout(2,2);
 
    nexttile;
    plot(time_signal, z_float_all{m}, 'k', 'LineWidth', 1);
    title(sprintf('%s Floating-point (64-bit)', method_names{m}));

    ylim([-4 4]); xlabel('Time (s)'); ylabel('Z-score'); grid on;
 
    for k = 1:total_num_widths
        nexttile;
        plot(time_signal, z_float_all{m}, 'k', 'LineWidth', 0.5); hold on;
        plot(time_signal, z_fixed_all{m}{k}, colors{k}, 'LineWidth', 1);
        hold off;
        title(sprintf('%s %d-bit vs. float',method_names{m}, total_word_lengths(k)));
        ylim([-4 4]); xlabel('Time (s)'); ylabel('Z-score'); grid on;
    end
end
 
%% Figure 5: Quantization error, small multiples -- method x word length
% Dense but complete: every combination gets its own tightly-scaled axes.
error_all = {error_fixed_direct, error_fixed_mitchell, error_fixed_pla};
 
figure('Name', 'Quantization error, all methods x word lengths');
tiledlayout(3, total_num_widths);
 
for m = 1:3
    for k = 1:total_num_widths
        nexttile;
        plot(time_signal, error_all{m}{k}, colors{k}, 'LineWidth', 0.5);
        title(sprintf('%s, %d-bit', method_names{m}, total_word_lengths(k)));
        xlabel('Time (s)'); ylabel('Error');
        grid on;
    end
end
 
%% Figure 6: Quantization error vs. word length, all methods overlaid
figure('Name', 'Quantization Error Comparison Across Methods');
semilogy(total_word_lengths, mse_fixed_direct, '-o', 'LineWidth', 2, 'Color', 'r'); hold on;
semilogy(total_word_lengths, mse_fixed_mitchell, '-s', 'LineWidth', 2, 'Color', 'b');
semilogy(total_word_lengths, mse_fixed_pla, '-^', 'LineWidth', 2, 'Color', 'g');
hold off;
xlabel('Fixed-point word length (bits)');
ylabel('MSE (log scale)');
title('Fixed-point quantization error by method');
legend('Direct', 'Mitchell', 'PLA', 'Location', 'best');
grid on;
 
%% Figure 7: Floating-point approximation accuracy (vs. direct as reference)
accuracy_gap = [0, ...   % Direct compared to itself: trivially zero, included for symmetry
                mean(abs(z_mitchell_float - z_direct_float)), ...
                mean(abs(z_pla_float - z_direct_float))];
 
figure('Name', 'Floating-point approximation accuracy');
bar(categorical({'Direct (reference)', 'Mitchell', 'PLA'}), accuracy_gap);
ylabel('Mean |error| vs. direct-window ground truth (V)');
title('Inherent approximation error, independent of quantization');
grid on;
 
%% Figure 8: MSE grouped bar chart, method x word length, with value labels
mse_matrix = [mse_fixed_direct'; mse_fixed_mitchell'; mse_fixed_pla'];
 
figure('Name', 'MSE by method and word length');
b = bar(categorical({'8-bit','16-bit','32-bit'}), mse_matrix', 'grouped');
set(gca, 'YScale', 'log');
ylabel('MSE (log scale)');
legend('Direct','Mitchell','PLA', 'Location', 'best');
title('Quantization error grouped by word length');
grid on;
 
for i = 1:numel(b)
    xtips = b(i).XEndPoints;
    ytips = b(i).YEndPoints;
    labels = arrayfun(@(v) sprintf('%.2g', v), b(i).YData, 'UniformOutput', false);
    text(xtips, ytips, labels, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 8, 'FontWeight', 'bold');
end

figure('Name', 'Method agreement');
plot(time_signal, z_mitchell_float - z_direct_float, 'b', 'LineWidth', 0.5); hold on;
plot(time_signal, z_pla_float - z_direct_float, 'g', 'LineWidth', 0.5);
hold off;
xlabel('Time (s)'); ylabel('Difference from Direct (Z-score units)');
legend('Mitchell - Direct', 'PLA - Direct');
title('Inter-method agreement over time');
grid on;

%% 8. Plotting Signal vs Power

% % figure('Name', 'Raw Indy');
% subplot(1,2,1);
% 
% plot(time_signal, voltage_signal, 'k', 'LineWidth', 2);
% %hold on;
% %plot(time_signal, MAT_teFIR, 'r', 'LineWidth', 2);
% %hold off;
% xlabel('Time (s)'); ylabel('Voltage (V)');
% %legend('Before (Unfiltered)', 'After (MATLAB FIR)');
% set(gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');
% grid on;
% 
% subplot(1,2,2);
% imagesc(t, f, s_dB);
% axis xy; 
% caxis(color_limits); 
% xlabel('Time (s)');
% ylabel('Frequency (Hz)'); 
% set(gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');
% 
% % --- Single Colorbar Configuration ---
% cb = colorbar; 
% cb.Label.String = 'Power Magnitude (dB)'; 
% cb.Label.FontSize = 18; 
% cb.Label.FontWeight = 'bold';


%% 9. Noise Generated Models

fs_lfp = 1000; % Sampling frequency (samples per second)
f_syn = 20; % Frequency of signal in Hz
t = linspace(0, 1, 1000); % 1000 ms of LFP data
a_peak = 200 * 10^((-6)); % Peak Amplitude
y = a_peak * sin(2 * pi * f_syn * t);

% % *White Noise (Thermal)
% fprintf('\n\nGenerating White Noise...\n');
% n_power = 0.0000000001; % Noise power
% noise_start = 1;
% noise_end = 1000;
% noise_length = noise_end - noise_start + 1;
% segment_white_noise = sqrt(n_power) * randn(numChannels,noise_length);
% signal_with_white = y + segment_white_noise;
% 
% % *OU Noise (Colored)
% fprintf('Generating OU Noise...\n');
% fs = 1000; % Sampling frequency in Hz
% tau_m = 0.020; % Membrane time constant: 20ms standard (Pyramidal cells)
% theta = 1 / tau_m; % Reversion rate (cutoff frequency)
% u = 0; % Resting potential
% sigma = 1; % Noise amplitude (intensity)
% dt = 1 / fs; % Time step
% ou_noise = zeros(size(y)); 
% ou_noise(1) = randn * sigma;
% 
% for i = 1:length(y)-1
%     drift = theta * (u - ou_noise(i)) * dt;
%     diffusion = sigma * sqrt(dt) * randn;   % randn (xi)
%     ou_noise(i+1) = ou_noise(i) + drift + diffusion;
% end
% ou_noise = 0.0001 * ou_noise;
% signal_with_ou_noise = y + ou_noise; % Synthetic Noise

% % *Motion artifacts
% fprintf('Generating Motion Noise...\n');
% index_shift = 0.5 * fs_lfp;
% decay_rate = -2;
% motion = zeros(size(t));
% motion(index_shift:end) = 1.5 * exp(decay_rate * (0:(length(t)-index_shift))/fs_lfp);
% motion = 0.0001 * motion; % Scale the amplitude
% signal_with_motion_noise = y + motion;
% 
% figure('Name', 'Raw Indy');
% subplot(1,2,1);
% 
% plot(t, y, 'k', 'LineWidth', 2);
% xlabel('Time (s)'); ylabel('Voltage (V)');
% %legend('Before (Unfiltered)', 'After (MATLAB FIR)');
% set(gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');
% grid on;
% ylim([-0.0003 0.0003]) 
% 
% subplot(1,2,2);
% plot(t, signal_with_motion_noise, 'k', 'LineWidth', 2);
% xlabel('Time (s)'); ylabel('Voltage (V)');
% %legend('Before (Unfiltered)', 'After (MATLAB FIR)');
% set(gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');
% grid on;
% ylim([-0.0003 0.0003]) 
