% SNR Proxy module
function [band_status, bandpower_ratio] = SNR_Proxy(data, fs, signal_band, noise_band)
     %{
        Description: Estimate the bandpower ratio of neural activity
        by relative strength vs background noise using a STFT.
       
        P_i(B) = sum_{f in B} |X_i(f)|^2
        SNR_i = P_i(B_sig) / (P_i(B_noise) + epsilon)

        Inputs:
            data - The 1D data vector
            fs - Sampling frequency (Hz)
            signal_band - Task relevant band [f_low f_high]
            noise_band - Noise relevant band [f_low f_high]

        Output:
            bandpower_ratio - The calculated ratio per-frame
            band_status - The channel conditional status (mean SNR check)
     %}  

    % Adding epsilon to avoid division by zero if noise is 0
    epsilon = 1e-6;

    % SNR threshold (linear ratio) for channel inclusion
    tau_snr = 1.2; % Signal is 20% higher than noise

    % Short-Time Fourier Transforms: X_i(f) per frame
    [s, f, ~] = stft (data, fs);
    p_spectrogram_v = abs(s(:,:,1)).^2; % |X(f)|^2
    
    % Frequency-bin for each band
    signal_p_indices_v = f >= signal_band(1) & f <= signal_band(2);
    noise_p_indices_v  = f >= noise_band(1)  & f <= noise_band(2);

    % P_i(B) for sum over freq bins per frame
    power_signal = sum(p_spectrogram_v(signal_p_indices_v, :), 1);
    power_noise = sum(p_spectrogram_v(noise_p_indices_v, :), 1);
 
    % SNR_i = P_i(B_sig) / (P_i(B_noise) + epsilon)
    bandpower_ratio = power_signal ./ (power_noise + epsilon);
    
   if mean(bandpower_ratio) > tau_snr 
        band_status = 'Included!';
   else    
        band_status = 'Excluded!';
   end
end