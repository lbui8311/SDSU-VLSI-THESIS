% Calculate Running Variance Function
function sig_var = calculate_variance(data, alpha)
    %{
        Description: Calculating the EMA-based running variance
        Inputs:
            data - The 1D data vector
            alpha - The fixed alpha constant
        Output:
            sig_avg - The 1D calculated variance vector

        Equations:
            u(t) = a * x(t) + (1-a) * u(t-1)
            v(t) = a * (x(t) - u(t))^2 + (1-a) * v(t-1)
    %}

    % Variables declaration
    num_samples = length(data);
    u = zeros(num_samples, 1);
    x = data(:); % ensure vector

    % First value initialization
    u(1) = data(1);
    v(1) = 0;

    % Calculate mean and variance
    for t = 2:num_samples
        u(t) = alpha * x(t) + (1 - alpha) * u(t - 1);
        v(t) = alpha * (x(t) - u(t))^2 + (1 - alpha) * v(t - 1);
    end

    if isvector(v)
        v = v(:).'; % Transpose
    end

    sig_var = v(end); % Only grab newest value
    %sig_var = v;
end