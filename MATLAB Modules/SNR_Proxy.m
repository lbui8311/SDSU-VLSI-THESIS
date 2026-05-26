% SNR Proxy module
function [bandpower_ratio, band_status] = calculate_SNR(data)
     %{
        Description: Calculating the Signal-to-Noise Ratio
        Inputs:
            data - The 1D data vector
        Output:
            bandpower_ratio - The calcualted ratio
            band_status - The channel conditional status
     %}

    % Smoothing factor for running mean and variance
    alpha = 0.1;     
    % Adding epsilon to avoid division by zero if noise is 0
    epsilon = 1e-6;
    flag_list = []; % Keep track of all indices that fail
    
    % Calculating Noise Band using Median Absolute Deviation (MAD)
    noise_std_est = median(abs(data(1,:)))./0.6745; % Global noise 75th
    power_noise_est = noise_std_est.^2;
    
    % Calculating Signal Band using EMA-based Running Variance
    power_total = calculate_variance(data(1,:), alpha);
    power_signal =  max(power_total - power_noise_est, epsilon); % Ensure power signal is never negative

    bandpower_ratio = 10 * log10(power_signal(1,:) ./   (power_noise_est(1,:) + epsilon));
   
   if bandpower_ratio > 1.2 % Signal is 20% higher than noise
        band_status = 'Included!';
   else    
        band_status = 'Excluded!';
   end
end