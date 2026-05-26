% Hysteretic Masking Module
function [channel_status, channel_weight] = calculate_hysteretic(data)
     %{
        Description: Calculating the channel weights
        Inputs:
            data - The 1D data vector
        Output:
            channel_weight - The calcualted channel weights
            channel_status - The channel conditional status
     %}

    alpha = 0.1; % Smoothing factor for running mean and variance  
    W = 10;
    weight_list = [];
    total_fail_count = 0;
    total_pass_count = 0;
    channel_weight = 0;
    max_fail = 5;
    max_pass = 10;
    lambda = 0.9; % Weight decay factor if fail
    delta = 0.05; % Weight recover factor if pass
    tau_min = 1e-13; % Sample Amplifier noise floor ~0.3 uV RMS (Dead/Shorted)
    tau_max = 1e-7; % Intracortical LFP rarely exceeds ~300 uV RMS (Massive Artifacts)
    
    for i = 1:length(data)
        % Create a trailing window
        begin_index = max(1, i - W + 1); % max(1,...) for edge cases if negative
        end_index = i; % if i = 20, begin_index = 11, end_index = 20
          
        % Find mean and variance of current frame
        frame = data(1,begin_index:end_index);
        current_frame = frame(1,:);
        mean = calculate_mean(current_frame, alpha);
        variance = calculate_variance(current_frame, alpha); % in V^2
        
        if (variance > tau_max || variance < tau_min)
            total_fail_count = total_fail_count + 1;
            total_pass_count = 0;
        else
            total_pass_count = total_pass_count + 1;
            total_fail_count = 0;
        end

        % External Frame Limits
        if total_fail_count >= max_fail % Decaying factor
            channel_weight = channel_weight * lambda;
            total_fail_count = 0; 
        elseif total_pass_count >= max_pass % Recovery incrementing
            channel_weight = min(1, channel_weight + delta); 
            total_pass_count = 0; 
        else 
            channel_weight = channel_weight;
        end
        weight_list = [weight_list, channel_weight];
    end
    
    if channel_weight == 1
        channel_status = 'Healthy!';
    elseif channel_weight == 0
        channel_status = 'Unhealthy!';
    else
        channel_status = 'Transitioning!';
    end

end