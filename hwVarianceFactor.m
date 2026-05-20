function f = hwVarianceFactor(a, tau)
%% HWVARIANCEFACTOR  Closed-form exponential variance integral.
%
%  f(a,tau) = (1-exp(-2*a*tau))/(2*a), with limit tau as a -> 0.

if abs(a) < 1e-10
    f = tau;
else
    f = (1 - exp(-2 * a * tau)) / (2 * a);
end
end
