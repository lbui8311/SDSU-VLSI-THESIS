% Amplitude Clipping Module
function [channel_status, flag_list] = amplitude_clip(data)
    %{
        Description: Keep track of all flagged data point's amplitudes
        Inputs:
            data - The 1D data vector
        Output:
            channel_status - The channel conditional status
            flag_list - A list of all data points that got flagged
    %}

    W = 10; % Num of data points in a sliding window
    flag_idx = 0; % Flag index of exceeded amplitude
    flag_list = []; % Keep track of all flagged indices
    flag_limit = 100; % Num of permitted flagged indices
    max_amp = 100 * 10^(-6); % Maximum Amplitude point
    channel_status = 'Unchecked!';
    
    for i = 1:length(data)
        % Create a trailing window
        begin_index = max(1, i - W + 1); % max(1,...) for edge cases if negative
        end_index = i; % if i = 20, begin_index = 11, end_index = 20
        window = data(1,begin_index:end_index); % Grab current window

        % Finding peak amplitude(s) of current window
        for j = 1:length(window)
           if window(j) > max_amp
               flag_idx = i;
               flag_list = [flag_list, flag_idx];
           end
        end
    end
    
    % Final Channel Verdict
    if length(flag_list) > flag_limit
        channel_status = 'Unreliable!';
    else
        channel_status = 'Reliable!';
    end
end