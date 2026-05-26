function P = bondDirtyModel(bond, OIS, ZS, t0)
% BONDDIRTYMODEL  Dirty model price of a single corporate coupon bond.
%
% Input:
%       bond - struct for a SINGLE bond with fields
%              .couponDates : column of future payment dates (datenums)
%              .cashflows   : column of cash-flows (per unit face;
%                             last entry already includes the +1 redemption,
%                             as produced by loadBondData)
%       OIS  - risk-free discount curve with fields
%              .T  : column of pillar dates (datenums)
%              .DF : column of discount factors P^D(t0, T)
%       ZS   - Z-spread curve from bootstrapZSpread with fields
%              .T  : column of pillar dates (= bond maturities)
%              .Z  : column of Z-spreads (decimal, cont. comp., ACT/365)
%       t0   - settlement date (datenum)
%
% Output:
%       P    - dirty model price per unit face value

t_cf = bond.couponDates;
c_cf = bond.cashflows;

tau = yearfrac(t0, t_cf, 3);   % ACT/365
dfOIS = interpolateOIS(t0, OIS.T, OIS.DF, t_cf);
Z_cf = interpolateZSpread(t0, ZS.T, ZS.Z, t_cf);

Bbar = dfOIS .* exp(-Z_cf .* tau);
P = sum(c_cf .* Bbar);

end
