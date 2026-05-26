function z_score = calculate_z_score(data, alpha, epsilon)
    %{
        Description: Calculating the EMA-based Z-score normalization
        Inputs:
            data - The 1D data vector
            alpha - The fixed alpha constant
            epsilon - The fixed epsilon constant
        Output:
            z_score - The 1D calculated z-score vector

        Equations:
            u(t) = a * x(t) + (1-a) * u(t-1)
            v(t) = a * (x(t) - u(t))^2 + (1-a) * v(t-1)
            z(t) = (x(t) - u(t))/ sqrt(v(t)+e)
    %}

    % Variables declaration
    num_samples = length(data);
    u = zeros(num_samples, 1);
    v = zeros(num_samples, 1);

    x = data(:); % ensure vector
    
    % First values initialization
    u(1) = data(1);
    v(1) = 0;

    % Calculate mean and variance
    for t = 2:num_samples
        u(t) = alpha * x(t) + (1 - alpha) * u(t - 1);
        v(t) = alpha * (x(t) - u(t))^2 + (1 - alpha) * v(t - 1);
    end

    % Calculate z-score
    z = (x - u) ./ sqrt(v + epsilon);

    if isvector(z)
        z = z(:).'; % Transpose
    end

    z_score = z;
end