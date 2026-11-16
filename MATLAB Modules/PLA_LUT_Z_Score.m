% PLA LUT Z-score Fixed-point 6/8/16-bit
function [z_result_PLA, segment_idx_trace] = PLA_LUT_Z_Score(data, alpha, epsilon)
     %{
        Description: Calculating the PLA Z-score normalization, a
        piecewise-linear approximation of 1/sqrt(v + epsilon).

        Equations:

            EMA-based Mean and Variance:
            u(t) = alpha * x(t) + (1 - alpha) * u(t - 1);
            v(t) = alpha * (x(t) - u(t))^2 + (1 - alpha) * v(t - 1);

            norm: x(t) - u(t)
            denom: a_value * (v(t) + epsilon) + b_value or y = ax + b
                    
        Inputs:
            data - The 1D data vector
            W - The window size
            alpha - The fixed alpha constant
            epsilon - The fixed epsilon constant

        Output:
            z_result_PLA - The 1D calculated z-score vector

     %}

    numSamples = length(data);
    x = data(:).'; % Ensure vector

    z_result_PLA = zeros(1, length(data));

    % Creating vector x with number of points
    % xmin = epsilon; % Can't be 0 => inf
    % xmax = 0.0215; % max of var + eps
    xmin = 1e-12; % min var val float
    xmax = 5e-8; % max var val float
    %xLength = xmax - xmin;
    numSegments = 10;
    %segmentSize = xLength/numSegments;
    
    numPoints_X = length(data);
    %x_space = linspace(xmin, xmax, numPoints_X);
    x_space = logspace(log10(xmin), log10(xmax), numPoints_X);

    segmentCoeff = zeros(2, numSegments);
    xSplit = mat2cell(x_space, 1, ones(1, numSegments) .* length(x_space)/numSegments)';
    
    % Calculate and divide segments for 1 / sqrt(x)
    for j = 1 : numSegments
        cX = xSplit{j};
        cY = 1./sqrt(cX);
        p = polyfit(cX, cY, 1);
        segmentCoeff(:, j) = p';   
    end
    
    % Store a and b coefficients
    LUT_a = segmentCoeff(1,:);
    LUT_b = segmentCoeff(2,:);
    %points_per_segment = length(data) / numSegments;

    % var_float_trace = zeros(1, length(data));   % preallocate before the loop
    % segment_idx_trace = zeros(1, length(data));  % Test 

    u = zeros(1, numSamples);
    v = zeros(1, numSamples);

    % First values initialization
    u(1) = x(1);
    v(1) = 0;

    segment_idx_trace = zeros(1, numSamples);

    % Looping through all data points to find the current frame
    for t = 1:numSamples
        % % Create a trailing window
        % begin_index = max(1, i - W + 1); % max(1,...) for edge cases if negative
        % end_index = i; % if i = 20, begin_index = 11, end_index = 20
        % 
        % % Find mean and variance of current frame
        % frame = data(1,begin_index:end_index);
        % current_frame = frame(1,:);
        % 
        % mean_float = EMA_Running_Mean(current_frame, alpha); % Double 64-bit
        % var_float = EMA_Running_Variance(current_frame, alpha); % in V^2

        %var_float_trace(t) = var_float;

        % Continuous signal instead of window
        if t > 1
            u(t) = alpha * x(t) + (1 - alpha) * u(t - 1);
            v(t) = alpha * (x(t) - u(t))^2 + (1 - alpha) * v(t - 1);
        end

        mean_temp = u(t); % in V
        var_temp = v(t); % in V^2
    
        % Find current segment id
        val = var_temp + epsilon;
        %segment_idx = ceil(val / points_per_segment);
        %segment_idx = ceil((val - xmin) / segmentSize);
        segment_idx = ceil((log10(val) - log10(xmin)) / (log10(xmax) - log10(xmin)) * numSegments);
        segment_idx = max(1, min(segment_idx, numSegments)); % Prevent out of bounds
        segment_idx_trace(t) = segment_idx;   % Test

        % Locate the current a and b coefficients
        a_value = LUT_a(segment_idx);
        b_value = LUT_b(segment_idx);
    
        nom = double(x(t) - mean_temp);
        denom = a_value * (var_temp + epsilon) + b_value;
        z_result_PLA(t) = nom .* denom;
    end
end

