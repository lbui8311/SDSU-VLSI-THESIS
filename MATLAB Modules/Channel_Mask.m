function [masked_x_n, m_i_all, m_i_am, m_i_var, channel_weights, var_flags_sample] = Channel_Mask(x_n, tau_sat, var_min, var_max)

    %{
        Description: Evaluating the recording channel using Amplitude
        Saturation, Variance Stability and Hysteresis with thresholds.
       
        Equations:

            Amplitude Saturation Detection:
            A_i = |x_i[n]|
            S_i = {1,  A_i <= tau_sat
                  {0,  A_i > tau_sat
            
            Variance-Stability Detection (Normal)
            u_i = (1/N_M) * SUM_n=1_(N_M) x_i[n]
            o^2_i = (1/N_M) * SUM_n=1_(N_M) (x_i[n] - u_i)^2
            v_i = {1,  o^2_min <= o^2_i <= o^2_max
                  {0,  otherwise

            Variance-Stability Detection (Hardware Efficient)
            S_i = SUM_n=1_(N_M) x_i[n]
            Q_i = SUM_n=1_(N_M) x^2_i[n]
            u_i = S_i / N_M

            o^2_i = (1/N_M) * SUM_n=1_(N_M) (x^2_i[n] - u^2_i)
            or
            o^2_i = (N_M * Q_i - S^2_i) / N^2_M

            Hysteresis with thresholds:
            w_i = {lambda * w_i,           if fail for N consecutive frames,
                  {min(1, w_i + delta,     if pass for N consecutive frames

        Inputs:
            x_n - The 2D recording channels and samples
            tau_sat - The predefined saturation threshold
            var_min - The lower variance threshold
            var_max - The upper variance threshold

        Output:
            masked_x_n - The masked 2D recording channels and samples
            m_i_am - The channel validity status for Amplitude Saturation
            m_i_var - The channel validity status for Variance Stability
            m_i_all - The channel validity status for both am & var
    %}


    [num_channels, num_samples] = size(x_n);

    % Amplitude-Saturation F_sat variables for samples and channels
    s_i_am_samples = zeros(num_channels, num_samples);
    s_i_am_channels = zeros(num_channels, 1);
    
    % First test: Amplitude-Saturation Detection
    % Continuous monitoring
    for c = 1:num_channels
        for s = 1:num_samples
            if abs(x_n(c,s)) > tau_sat
                s_i_am_samples(c,s) = 0;
            else
                s_i_am_samples(c,s) = 1; % Good
            end
        end
        
        % If at least one sample is greater than threshold => 0 or invalid
        if max(abs(s_i_am_samples(c,:))) == 0
            s_i_am_channels(c) = 0;
        else
            s_i_am_channels(c) = 1; % Good
        end
    end
    
    % Final status of all channels
    m_i_am = s_i_am_channels;

    % Second test: Variance-Stability Detection
    % Normal approach setup variables
    temp_sum_mean = zeros(num_channels, num_samples);
    temp_sum_var = zeros(num_channels, num_samples);
    channels_mean = zeros(num_channels, num_samples);
    channels_var = zeros(num_channels, num_samples);

    var_flags_sample = zeros(num_channels, num_samples);
    var_flags_channel = zeros(num_channels, 1);

    % Hardware-efficient approach setup variables (optional)
    s_i_var = zeros(num_channels, num_samples);
    q_i_var = zeros(num_channels, num_samples);
    u_i_var = zeros(num_channels, num_samples);

    % Hysteresis variables (optional)
    max_fail = 5; % Consecutive fails before decay
    max_pass = 10; % Consecutive pass before recovery
    lambda = 0.9; % Weight decay factor if fail
    delta = 0.05; % Weight recover factor if pass
    channel_weights = zeros(num_channels,num_samples);
    total_fail_count = zeros(num_channels,1);
    total_pass_count = zeros(num_channels,1);

    % 
    % for c = 1:num_channels
    %     temp_weights = 0; total_pass_count = 0; total_fail_count = 0;
    %     for s = 1:num_samples
    % 
    %         % % Normal
    %         % if s == 1; temp_sum_mean(c,s) =  x_n(c,s);
    %         % else; temp_sum_mean(c,s) = x_n(c,s) + x_n(c,s-1); end
    %         % channels_mean(c,s) = (1/length(x_n(1,:)) * temp_sum_mean(c,s));
    %         % 
    %         % if s == 1; temp_sum_var(c,s) =  x_n(c,s);
    %         % else; temp_sum_var(c,s) = power((x_n(c,s) - channels_mean(c,s)), 2) + power((x_n(c,s-1) - channels_mean(c,s-1)), 2) ; end
    %         % channels_var(c,s) = (1/length(x_n(1,:)) * temp_sum_var(c,s));  
    % 
    %         % % Hardware-efficient (optional)
    %         if s == 1
    %             s_i_var(c,s) = x_n(c,s); 
    %             q_i_var(c,s) = power(x_n(c,s), 2);
    %         else 
    %             s_i_var(c,s) = x_n(c,s) + s_i_var(c,s-1);
    %             q_i_var(c,s) = power(x_n(c,s), 2) + q_i_var(c,s-1);
    %         end
    % 
    %         u_i_var(c,s) = s_i_var(c,s) ./ length(x_n(1,:));
    %         %channels_var(c,s) = (q_i(c,s) / length(x_n(1,:))) - power((s_i(c,s) / length(x_n(1,:))),2);
    %         channels_var(c,s) = (length(x_n(1,:)) * q_i_var(c,s)  - power(s_i_var(c,s),2)) ./  power(length(x_n(1,:)), 2);
    % 
    %         % Hysteretic Masking (optional)
    %         if (channels_var(c,s) > var_max || channels_var(c,s) < var_min)
    %             total_fail_count = total_fail_count + 1;
    %             total_pass_count = 0;
    %         else
    %             total_pass_count = total_pass_count + 1;
    %             total_fail_count = 0;
    %         end
    % 
    %         % External Frame Limits
    %         if total_fail_count >= max_fail % Decaying factor
    %             temp_weights = temp_weights * lambda;
    %             total_fail_count = 0; 
    %         elseif total_pass_count >= max_pass % Recovery increment
    %             temp_weights = min(1, temp_weights + delta); 
    %             total_pass_count = 0; 
    %         end
    % 
    %         channel_weights(c,s) = temp_weights;
    % 
    %         % %if channels_var(c,s) >= var_min && channels_var(c,s) <= var_max
    %         % if channel_weights(c) >= var_min && channel_weights(c) <= var_max
    %         %     var_flags_sample(c,s) = 1; % Good
    %         % else
    %         %     var_flags_sample(c,s) = 0;
    %         % end
    %     end
    % 
    %     %if channels_var(c,s) >= var_min && channels_var(c,s) <= var_max
    %     var_flags_sample = channel_weights == 1; % 1 is good
    % 
    %     % If the sample weight has never fully recovered then channel is invalid
    %     var_flags_channel(c) = any(var_flags_sample(c,:)); % 1 is good
    % end

    % s_i_var = sum(x_n(:, 1:num_samples), 2);
    % q_i_var = sum(x_n(:, 1:num_samples).^2, 2);
    % %u_i_var = s_i_var ./ num_samples;
    % channels_var = (q_i_var / num_samples) - (s_i_var / num_samples).^2;

    % Variables declaration
    u = zeros(num_channels, num_samples);
    v = zeros(num_channels, num_samples);

    % First value initialization
    u(:,1) = x_n(:,1);
    v(:,1) = 0;

    % Calculate EMA-based mean and variance
    for c = 1:num_channels
        for t = 2:num_samples
            u(c,t) = 0.1 * x_n(c,t) + (1 - 0.1) * u(c,t - 1);
            v(c,t) = 0.1 * (x_n(c,t) - u(c,t))^2 + (1 - 0.1) * v(c,t - 1);
        end
    end

    bad_flags = v > var_max | v < var_min; % 10 x 100000 logical: 1 = this sample fails (bad)
    good_flags = ~bad_flags;  % same, but 1 = good 

    var_flags_sample = bad_flags;
    %disp(size(var_flags_sample));
    
    %channels_var = sum(x_n(:, 1:num_samples).^2 - u_i_var(:, 1:num_samples).^2 , 2) ./ num_samples;

    %channels_var = ((num_samples) .* q_i_var  - power(s_i_var,2)) ./  power(num_samples, 2);

    m_i_var = var_flags_channel;

    % Combine amplitude and variance checks to create the channel mask
    m_i_all = m_i_am .* m_i_var;
    masked_x_n = x_n .* m_i_all;
end