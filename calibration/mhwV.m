function v = mhwV(a, zeta_alpha, tau)
%% MHWV  Coefficient v_{alpha,i}.
%
%  v_{alpha,i} = zeta_alpha * (1 - exp(-a*tau)) / a,
%  where tau = t_i - t_alpha.

if abs(a) < 1e-10
    v = zeta_alpha * tau;
else
    v = zeta_alpha * (1 - exp(-a * tau)) / a;
end
end
