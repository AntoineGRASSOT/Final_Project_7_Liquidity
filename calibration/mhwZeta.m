function zeta = mhwZeta(a, sigma, t_alpha)
%% MHWZETA  Standard deviation of the MHW Gaussian factor.

if abs(a) < 1e-10
    zeta2 = sigma^2 * t_alpha;
else
    zeta2 = sigma^2 * (1 - exp(-2 * a * t_alpha)) / (2 * a);
end

zeta = sqrt(zeta2);
end
