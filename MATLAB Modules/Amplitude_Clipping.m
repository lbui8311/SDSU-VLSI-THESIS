% Amplitude Clipping Module
function capped_signal = Amplitude_Clipping(data, max_amp)
    %{
        Description: Take in the filtered signal and set all data points that 
        exceeded plausible voltage range to maximum amplitude.

        Inputs:
            data - The 1D data vector
            max_amp - The maximum amplitude or out of range limit
            
        Output:            
            capped_signal - The newly capped to max amplitude signal
    %}

    W = 10; % Num of data points in a sliding window
    [num_channels, num_samples] = size(data); 
    capped_signal = zeros(num_channels, num_samples);
    
    for i = 1:num_samples
        % Create a trailing window
        begin_index = max(1, i - W + 1); % max(1,...) for edge cases if negative
        end_index = i; % if i = 20, begin_index = 11, end_index = 20
        window = data(1,begin_index:end_index); % Grab current window

        % Finding peak amplitude(s) of current window
        if max(abs(window)) > max_amp
            capped_signal(1,i) = sign(data(i)) * max_amp;
        else 
            capped_signal(1,i) = data(i);
        end
    end
    
end