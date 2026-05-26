% PLA LUT Z-score Fixed-point 6/8/16-bit
function z_result_PLA = z_score_PLA(data, W, alpha, epsilon)
     %{
        Description: Calculating the PLA Z-score normalization
        Inputs:
            data - The 1D data vector
            W - The window size
            alpha - The fixed alpha constant
            epsilon - The fixed epsilon constant
        Output:
            z_result_PLA - The 1D calculated z-score vector
        Equations:
            norm = x(t) - u(t)
            denom = a_value * (v(t) + epsilon) + b_value % y = ax + b
     %}

    z_result_PLA = zeros(1, length(data));

    % Creating vector x with number of points
    xmin = epsilon; % Can't be 0 => inf
    xmax = 0.0215; % max of var + eps
    xLength = xmax - xmin;
    numSegments = 10;
    segmentSize = xLength/numSegments;
    
    numPoints_X = length(data);
    x = linspace(xmin, xmax, numPoints_X);
    
    segmentCoeff = zeros(2, numSegments);
    xSplit = mat2cell(x, 1, ones(1, numSegments) .* length(x)/numSegments)';
    
    % Calculate and divide segments for 1 / sqrt(x)
    for j = 1 : numSegments
        cX = xSplit{j};
        cY = 1./sqrt(cX);
        p = polyfit(cX, cY, 1);
        segmentCoeff(:, j) = p';   
    end
    
    % Store a and b coeffients
    LUT_a = segmentCoeff(1,:);
    LUT_b = segmentCoeff(2,:);
    points_per_segment = length(data) / numSegments;

    % Looping through all data points to find the current frame
    for i = 1:length(data)
        % Create a trailing window
        begin_index = max(1, i - W + 1); % max(1,...) for edge cases if negative
        end_index = i; % if i = 20, begin_index = 11, end_index = 20
    
        % Find mean and variance of current frame
        frame = data(1,begin_index:end_index);
        current_frame = frame(1,:);
    
        mean_float = calculate_mean(current_frame, alpha); % Double 64-bit
        var_float = calculate_variance(current_frame, alpha); % in V^2
    
        % Find current segment id
        val = var_float + epsilon;
        segment_idx = ceil(val / points_per_segment);
        segment_idx = max(1, min(segment_idx, numSegments)); % Prevent out of bounds
        
        % Locate the current a and b coefficients
        a_value = LUT_a(segment_idx);
        b_value = LUT_b(segment_idx);
    
        nom = double(data(i) - mean_float);
        denom = a_value * (var_float + epsilon) + b_value;
        z_result_PLA(i) = nom .* denom;
    end
end