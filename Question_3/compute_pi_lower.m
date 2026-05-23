function pi_low = compute_pi_lower(ti, tN, tau, sigma, a)
    % Compute individual ksi values
    ksi_i = compute_ksi(sigma, a, ti, tau);
    ksi_N = compute_ksi(sigma, a, tN, tau);
    
    % notation : Sigma_i(tau) and Sigma_N(tau) are the cumulative vol
    S_i = cumul_vol(ksi_i, a, tau);
    S_N = cumul_vol(ksi_N, a, tau); 
    
    integrand = @(n) (exp(-1/8 * S_N^2) ./ (pi * sqrt(1 - n) .* sqrt(n))) .* ...
                     exp(-n/2 * S_i * (S_i - S_N)) .* ...
                     (1 + sqrt(pi * (1 - n) / 2) * S_N * exp((1 - n)/8 * S_N^2) .* normcdf(sqrt((1 - n)/2) * S_N)) .* ...
                     (1 + sqrt(pi * n / 2) * (2*S_i - S_N) * exp(n/8 * (2*S_i - S_N)^2) .* normcdf(sqrt(n/2) * (2*S_i - S_N)));
                 
    % Avoid endpoint singularities with tight tolerances
    pi_low = integral(integrand, 0, 1, ...
        'AbsTol', 1e-10, 'RelTol', 1e-8, ...
        'Waypoints', [0.01, 0.5, 0.99]);
end