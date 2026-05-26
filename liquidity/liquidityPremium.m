function Delta = liquidityPremium(bond, OIS, ZS, MHWparams, t0, tau, bound)
% LIQUIDITYPREMIUM  Absolute sheer liquidity premium Delta_tau for a single bond.
%
% Input:
%   bond       struct for a SINGLE bond with fields
%              .couponDates : column of future payment dates (datenums)
%              .cashflows   : column of cash-flows (face = 1, last incl. +1)
%   OIS        risk-free curve with fields .T (pillar dates), .DF
%   ZS         Z-spread curve with fields .T (pillar dates), .Z
%   MHWparams  vol params with .a and .sigma
%   t0         settlement date (datenum)
%   tau        time-to-liquidate date (datenum)
%   bound      'U' for upper bound (piUpper),
%              'L' for lower bound (piLower)
%
% Output:
%   Delta      scalar, sheer liquidity premium per unit face value

t_cf = bond.couponDates;
c_cf = bond.cashflows;
t_N = t_cf(end);    % Bond maturity (last cashflow date)  

% Keep only cashflows strictly after the liquidation date
idx = (t_cf > tau);
t_cf = t_cf(idx);
c_cf = c_cf(idx);

% Edge case: no cashflow after tau -> liquidity premium is zero by definition
if isempty(t_cf)
    Delta = 0;
    return;
end

[Sigma_i, Sigma_N] = cumulatedVol(MHWparams, tau, t_cf, t_N, t0);
dfOIS_cf = interpolateOIS(t0, OIS.T, OIS.DF, t_cf);
Z_cf = interpolateZSpread(t0, ZS.T, ZS.Z, t_cf);
tau_cf = yearfrac(t0, t_cf, 3);
B_cf = dfOIS_cf .* exp(-Z_cf .* tau_cf);

% Closed-form factor pi_i depending on the requested bound
if bound == 'U'
    pi_vec = piUpper(Sigma_i);
else
    pi_vec = arrayfun(@(s) piLower(s, Sigma_N), Sigma_i);
end

surv = survivalProb(ZS, t0, tau);

Delta = sum(c_cf .* B_cf .* (pi_vec - surv));

end
