function PD = getPD(EUR6M, T_query)
%% GETPD  Interpolate Euribor 6m pseudo-discount factors.
%
%  Interpolates log-pseudo-discount factors with pchip.

logPD = interp1(EUR6M.T, log(EUR6M.PD), T_query(:), 'pchip', 'extrap');
PD    = exp(logPD);

if isrow(T_query)
    PD = PD.';
end
end
