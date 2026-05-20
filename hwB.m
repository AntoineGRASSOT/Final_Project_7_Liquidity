function B = hwB(a, tau)
%% HWB  Hull-White coefficient B(a,tau) = (1-exp(-a*tau))/a.

if abs(a) < 1e-10
    B = tau;
else
    B = (1 - exp(-a * tau)) / a;
end
end
