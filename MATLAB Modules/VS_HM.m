% Hysteretic Masking Module
function [channel_status, channel_weight_list] = VS_HM(data, alpha, W, max_fail, max_pass, lambda, delta, tau_min, tau_max)
     %{
        Description: Calculating the per-channel weight using EMA-based
        mean and variance tracking with debouncing (Hysteretic Masking).

        w_i = {lambda * w_i,           if fail for N consecutive frames,
              {min(1, w_i + delta,     if pass for N consecutive frames

        Inputs:
            data - The 1D data vector
            alpha - Smoothing factor for running mean and variance 
            W - Sliding windows length in samples
            max_fail - Consecutive fails before decay
            max_pass - Consecutive pass before recovery
            lambda - Decaying weight factor on fail-streak
            delta - Recovering weight step on pass-streak
            tau_min - Minimum acceptable EMA variance
            tau_max - Maximum acceptable EMA variance

        Output:
            channel_weight_list - The calculated channel weights
            channel_status - The channel conditional status
     %}

    total_fail_count = 0;
    total_pass_count = 0;
    channel_weight = 0;

    [~, num_samples] = size(data);  
    channel_weight_list = zeros(1, num_samples);
    
    for i = 1:length(data)
        % Create a trailing window
        begin_index = max(1, i - W + 1); % max(1,...) for edge cases if negative
        end_index = i; % if i = 20, begin_index = 11, end_index = 20
          
        % Find mean and variance of current frame
        frame = data(1,begin_index:end_index);
        current_frame = frame(1,:);
        var_val = EMA_Running_Variance(current_frame, alpha); % in V^2
        
        if (var_val > tau_max || var_val < tau_min)
            total_fail_count = tr + 1;
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
        %else channel_weight = channel_weight;
        end
        channel_weight_list(i) = channel_weight;
    end
    
    if channel_weight == 1
        channel_status = 'Healthy!';
    elseif channel_weight == 0
        channel_status = 'Unhealthy!';
    else
        channel_status = 'Transitioning!';
    end

end
