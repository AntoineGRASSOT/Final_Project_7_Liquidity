function vol = bachelierCS_implVol(S0, K, T_exp, B_alpha, C_S0, price_mkt)
%% BACHELIERCS_IMPLVOL  Normal implied volatility from a CS swaption price.
%
%  Numerically inverts bachelierCS to find sigma such that:
%    bachelierCS(S0, K, T_exp, B_alpha, C_S0, sigma) = price_mkt
%
%  Uses bisection on [0, 5]. Prices are monotone increasing in volatility.
%
%  Inputs/outputs are the same as bachelierCS, with price_mkt replacing vol_norm.
%  Output vol is the normal implied volatility in decimal form.

priceFun = @(sig) bachelierCS(S0, K, T_exp, B_alpha, C_S0, sig);

vol_lo = 1e-6;
vol_hi = 5.0;
tol    = 1e-12;

if priceFun(vol_hi) < price_mkt 
    warning('bachelierCS_implVol: price above intrinsic value');
    vol = vol_hi; return;
end

if priceFun(vol_lo) > price_mkt
    warning('bachelierCS_implVol: price below intrinsic value');
    vol = vol_lo; return;
end

for iter = 1:200
    vol_mid = (vol_lo + vol_hi) / 2;
    if priceFun(vol_mid) < price_mkt
        vol_lo = vol_mid;
    else
        vol_hi = vol_mid;
    end
    if (vol_hi - vol_lo) < tol, break; end
end
vol = (vol_lo + vol_hi) / 2;
end
