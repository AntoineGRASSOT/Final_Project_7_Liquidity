function DF = getDF(OIS, T_query)
%% GETDF  Interpolate OIS discount factors.
%
%  Interpolates log-discount factors with pchip, consistently with the MHW
%  swaption library.

logDF = interp1(OIS.T, log(OIS.DF), T_query(:), 'pchip', 'extrap');
DF    = exp(logDF);

if isrow(T_query)
    DF = DF.';
end
end
