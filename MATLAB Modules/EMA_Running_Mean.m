% Calculate Running Mean Function
function sig_avg = EMA_Running_Mean(data, alpha)
    %{
        Description: Calculating the EMA-based running mean.

        Equation: u(t) = a * x(t) + (1-a) * u(t-1)

        Inputs:
            data - The 1D data vector
            alpha - The fixed alpha constant
        Output:
            sig_avg - The 1D calculated mean vector

    %}
    
    % Variables declaration
    num_samples = length(data);
    u = zeros(num_samples, 1);
    x = data(:); % ensure vector

    % First value initialization
    u(1) = data(1);

    % Calculate mean
    for t = 2:num_samples
        u(t) = alpha * x(t) + (1 - alpha) * u(t - 1);
    end

    if isvector(u)
        u = u(:).'; % Transpose
    end

    sig_avg = u(end); % Only grab newest value
    %sig_avg = u;
end