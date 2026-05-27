function price = bachelierCS(S0, K, T_exp, B_alpha, C_S0, vol_norm)
%% BACHELIERCS  Bachelier (Normal-Black) price for a CS receiver swaption.
%
%  Market formula:
%    R^{C,mkt} = B(t0,t_alpha) * C(S0) * { [K-S0]*N(-d) + sigma*sqrt(T)*phi(d) }
%    d = (S0-K)/(sigma*sqrt(T))
%
%  Inputs:
%    S0       - ATM swap rate, decimal
%    K        - strike rate
%    T_exp    - option expiry in years
%    B_alpha  - OIS discount factor B(t0,t_alpha)
%    C_S0     - cash annuity C_{alpha,omega}(S0)
%    vol_norm - normal volatility, decimal, e.g. 0.00647 = 64.7 bps
%
%  Output:
%    price - CS receiver swaption price per unit notional

sqrt_T = sqrt(max(T_exp, 1e-10));
sigma_tot = vol_norm * sqrt_T;

if sigma_tot < 1e-12
    price = B_alpha * C_S0 * max(K - S0, 0);
    return;
end

d     = (S0 - K) / sigma_tot;
price = B_alpha * C_S0 * ((K - S0)*normcdf(-d) + sigma_tot*normpdf(d));
end
