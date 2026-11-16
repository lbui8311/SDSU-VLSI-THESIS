% Mitchell's Approximation
function z_result_mitchell = Mitchell_Approx_Z_Score(data, alpha, epsilon)
    %{
        Description: Calculating the Mitchell Z-score normalization via log
        domain, logarithmic approximation of division.
       
        Equations:

            EMA-based Mean and Variance:
            u(t) = alpha * x(t) + (1 - alpha) * u(t - 1);
            v(t) = alpha * (x(t) - u(t))^2 + (1 - alpha) * v(t - 1);

            norm: x(t) - u(t)
            denom: (v(t) + epsilon) / 2

        Inputs:
            data - The 1D data vector
            W - The window size
            alpha - The fixed alpha constant
            epsilon - The fixed epsilon constant

        Output:
            z_result_mitchell - The 1D calculated z-score vector
    %}

    numSamples = length(data);
    x = data(:).'; % Ensure vector

    z_result_mitchell = zeros(1, numSamples);
    u = zeros(1, numSamples);
    v = zeros(1, numSamples);

    % First values initialization
    u(1) = x(1);
    v(1) = 0;
    
    % Looping through all data points to find the current frame
    for t = 1:numSamples

        % % Create a trailing window
        % begin_index = max(1, i - W + 1); % max(1,...) for edge cases if negative
        % end_index = i; % Ex: if i = 20, begin_index = 11, end_index = 20
        % 
        % % Find current frame
        % frame = data(1,begin_index:end_index);
        % current_frame = frame(1,:);

        % % Find mean and variance of current frame
        % mean_float = calculate_mean(current_frame, alpha); % Double 64-bit
        % var_float = calculate_variance(current_frame, alpha); % in V^2

        % Continuous signal instead of window
        if t > 1
            u(t) = alpha * x(t) + (1 - alpha) * u(t - 1);
            v(t) = alpha * (x(t) - u(t))^2 + (1 - alpha) * v(t - 1);
        end

        mean_temp = u(t); % in V
        var_temp = v(t); % in V^2

        % Perform nominator and denominator approximations
        nom = double(x(t) - mean_temp);
        denom = var_temp + epsilon;
        sign_p = sign(nom) .* sign(denom); % Extract the sign
        nom = max(abs(nom), epsilon); % to avoid -InF and NaN for 0
        if nom == 0
            nom = mean_float;
        end
        denom = max(abs(denom), epsilon); % Without change signal dynamics
    
        % Mitchel Approximation
        % Find k integer, use priority encoder in hardware to find MSB
        k_a = floor(log2(abs(nom)));  % Ex: log(7) ~ 2.8 = 2
        k_b = floor(log2(abs(denom)));  % Ex:  log(3) ~ 1.6 = 1
    
        % Find x fraction, shift left by (WL-1)-k or fast barrel shifter in HW
        % Take WFO bit as needed, use x = 2^E * M for binary number
        x_a = (abs(nom) ./ 2.^k_a)-1; % Ex: (7 / 2^2) - 1 = 1.75 - 1 = 0.75
        x_b = (abs(denom) ./ 2.^k_b)-1; % Ex: (3 / 2^1) - 1 = 1.5 - 1 = 0.5
    
        % Add/sub logarithms, calculate k and x of product and quotient
        log_a_approx = k_a + x_a; % Ex: 2.75
        log_b_approx = (k_b + x_b) / 2; % Ex:  1.5
        log_p_approx_div = log_a_approx - log_b_approx; % Ex:  1.25
    
        k_p_div = floor(log_p_approx_div); % Ex: 1
        x_p_div = log_p_approx_div - k_p_div; % Ex: 1.25 - 1 = 0.25
    
        % Get the magnitude of approx product & quotient
        z_result_mitchell(t) = sign_p .* (2.^k_p_div) .* (1 + x_p_div);
    end
end
