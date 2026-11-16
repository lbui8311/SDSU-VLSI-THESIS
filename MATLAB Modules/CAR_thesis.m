% Clear Workspace, commands and windows
clear; clc; close all;

%% 1. Reading BCI Data
% Define Data & Timestamps Paths
data_file = 'C:\Users\Dr_ Luan\Documents\MATLAB\E218-BCI-Files\Sabes_Data\dataset_files\indy_20160624_03_fs10k.nwb';
data_read_path = '/acquisition/timeseries/broadband/data';
timestamps_read_path = '/acquisition/timeseries/broadband/timestamps';

f  = 'C:\Users\Dr_ Luan\Documents\MATLAB\E218-BCI-Files\Sabes_Data\dataset_files\mearec_test_10s.h5';
X  = h5read(f, '/recordings').';                 % transpose to samples x channels
f_s = double(h5read(f, '/info/recordings/fs'));
P  = h5read(f, '/channel_positions').';
st0 = h5read(f, '/spiketrains/0/times');          % neuron 0; folders are numbered 0, 1, 2, ...

% Read the Noisy_channels & Power_Noise_Channels
noisy_read = h5read(data_file, '/info/noisy_channels');
power_read = h5read(data_file, '/info/power_noise_channels');

% Update data retention to 60 seconds (matching numSamples = fs * 60)
fs = 10e3;
numSamples = fs * 10; % 10 instead of 60 for quick test purpose
numChannels = 10;
filtered_data = zeros(numChannels, numSamples);

% Read the Data & Timestamps
voltage_signal = h5read(data_file, data_read_path);
voltage_signal = voltage_signal(1:numChannels,1:numSamples);
time_signal = h5read(data_file, timestamps_read_path);
time_signal = time_signal(:,1:numSamples);
 
for k = 1 : numChannels
    % Perform filtering step on all channels
    filtered_data(k, :) = FilterSBP(voltage_signal(k, :), 10e3, 300, 3000, 'filtfilt', 'order', 2);
end

%% 2. Masking the channels

% Channel Masking setup parameters
tau_sat = 1000 * 10^(-6); % Maximum Amplitude point
var_min = 1e-11; % Lower variance threshold
var_max = 1e-5; % Upper variance threshold
cv = 0;

% Running the channel-masking: Amplitude Saturation + Variance Stability
[masked_x_n, m_i_all, m_i_am, m_i_var, channel_weights] = Channel_Mask(filtered_data, tau_sat, var_min, var_max);

%% 3. Conventional CAR

data_mean_inst = mean(filtered_data, 1); % Find the inst avg across all channels
data_car_all = filtered_data - data_mean_inst; % Channel samples - inst avg

%% 3. Mask-aware CAR

% The number of currently valid channels:
for i = 1:length(m_i_all(:,1))
    if m_i_all(i) == 1
        cv = cv + 1;
    end
end

fprintf("\nNumber of valid channels: %d\n", cv);

% Normal mask-aware CAR variables initialization
s_v_temp = zeros(numChannels, numSamples); % for holding immediate sum values
x_bar_v = zeros(numChannels, numSamples);
x_car = zeros(numChannels, numSamples);

% Reciprocal mask-aware CAR variables initialization
r_LUT = zeros(numChannels, 1); % Reciprocal lookup table
x_bar_r = zeros(numChannels, numSamples);
x_car_r = zeros(numChannels, numSamples);

% Storing reciprocal coefficients to LUT 
for k = 1:numChannels
    r_LUT(k) = 1 / k;
end

% Calculating the reference using only valid-channel aggregate
for c = 1:numChannels
    for s = 1:numSamples
        if s == 1; s_v_temp(c,s) =  filtered_data(c,s) * m_i_all(c);
        else; s_v_temp(c,s) = (filtered_data(c,s) * m_i_all(c)) + s_v_temp(c,s-1); end
    end
    x_bar_v(c,:) = s_v_temp(c,:) ./ cv;
    x_bar_r(c,:) = s_v_temp(c,:) .* r_LUT(c);
end


% Calculate CAR output valid channels for normal and reciprocal
for c = 1:numChannels
    x_car(c, :) = m_i_all(c) .* (filtered_data(c, :) - x_bar_v(c, :));
    x_car_r(c, :) = m_i_all(c) .* (filtered_data(c, :) - x_bar_r(c, :));
end


%% 4. Plotting some channels

figure
tiledlayout(2,2);
nexttile;
plot(time_signal, filtered_data(6,:), 'm', 'LineWidth', 2);
hold on;
plot(time_signal, data_car_all(6,:), 'b', 'LineWidth', 2);
hold on;
plot(time_signal, x_car(6,:), 'r', 'LineWidth', 2);
hold on;
plot(time_signal, x_car_r(6,:), 'c', 'LineWidth', 2);
hold off;
xlabel('Time (s)'); ylabel('Voltage (V)');
title(sprintf('Channel %d out of %d', 6, numChannels));
legend('BP Filtered Data', 'Conv. CAR', 'Norm. M-A CAR', 'Recipr. M-A CAR');
set(gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');
grid on;

nexttile;
plot(time_signal, filtered_data(5,:), 'm', 'LineWidth', 2);
hold on;
plot(time_signal, data_car_all(5,:), 'b', 'LineWidth', 2);
hold on;
plot(time_signal, x_car(5,:), 'r', 'LineWidth', 2);
hold on;
plot(time_signal, x_car_r(5,:), 'c', 'LineWidth', 2);
hold off;
xlabel('Time (s)'); ylabel('Voltage (V)');
title(sprintf('Channel %d out of %d', 5, numChannels));
legend('BP Filtered Data', 'Conv. CAR', 'Norm. M-A CAR', 'Recipr. M-A CAR');
set(gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');
grid on;

nexttile;
plot(time_signal, filtered_data(4,:), 'm', 'LineWidth', 2);
hold on;
plot(time_signal, data_car_all(4,:), 'b', 'LineWidth', 2);
hold on;
plot(time_signal, x_car(4,:), 'r', 'LineWidth', 2);
hold on;
plot(time_signal, x_car_r(4,:), 'c', 'LineWidth', 2);
hold off;
xlabel('Time (s)'); ylabel('Voltage (V)');
title(sprintf('Channel %d out of %d', 4, numChannels));
legend('BP Filtered Data', 'Conv. CAR', 'Norm. M-A CAR', 'Recipr. M-A CAR');
set(gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');
grid on;

nexttile;
plot(time_signal, filtered_data(3,:), 'm', 'LineWidth', 2);
hold on;
plot(time_signal, data_car_all(3,:), 'b', 'LineWidth', 2);
hold on;
plot(time_signal, x_car(3,:), 'r', 'LineWidth', 2);
hold on;
plot(time_signal, x_car_r(3,:), 'c', 'LineWidth', 2);
hold off;
xlabel('Time (s)'); ylabel('Voltage (V)');
title(sprintf('Channel %d out of %d', 3, numChannels));
legend('BP Filtered Data', 'Conv. CAR', 'Norm. M-A CAR', 'Recipr. M-A CAR');
set(gca, 'FontSize', 18, 'FontAngle', 'italic', 'FontWeight', 'bold');
grid on;
