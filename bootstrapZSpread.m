function ZS = bootstrapZSpread(bonds, OIS, t0, name)
% BOOTSTRAPZSPREAD  Bootstrap the Z-spread curve Z(t0, Tk) from liquid corporate bonds 
%
% Input:
%   bonds  - per-issuer struct from loadBondData with fields
%            .maturity    : column of MATLAB datenums (sorted ascending)
%            .coupon      : column of annual coupon rates (decimal)
%            .cleanPrice  : column of clean prices (per unit face)
%            .accrual     : column of accruals (per unit face, ACT/ACT ICMA)
%            .couponDates : cell column of future payment-date vectors
%            .cashflows   : cell column of cash-flow vectors (face = 1)
%   OIS    - risk-free discount curve with fields
%            .T  : column of pillar dates
%            .DF : column of discount factors P^D(t0, T)
%   t0     - settlement date
%
% Output:
%   ZS.T     - column of pillar dates  (= bonds.maturity)
%   ZS.Z     - column of Z-spreads at pillars (decimal, cont. comp.)
%   ZS.DFbar - column of risky discount factors
%              Bbar(t0, T) = P^D(t0, T) * exp(-Z(T) * yearfrac(t0, T, 3))
%

mat = bonds.maturity;
cf = bonds.cashflows;
cpnDates = bonds.couponDates;
invoice = bonds.invoice;
N = length(mat);

Z = zeros(N, 1);
T = mat;

for k = 1:N
    c_cf   = cf{k};
    t_cf   = cpnDates{k};
    target = invoice(k);

    tau_cf   = yearfrac(t0, t_cf, 3);
    dfOIS_cf = interpolateOIS(t0, OIS.T, OIS.DF, t_cf);

    if k == 1
        % First bond: no previous pillars -> all cash-flows are "unknown"
        knownPV     = 0;
        idx_unknown = true(size(t_cf));
    else
        idx_known   = (t_cf <= T(k - 1));
        idx_unknown = ~idx_known;

        Z_known    = interpolateZSpread(t0, T(1:k-1), Z(1:k-1), t_cf(idx_known));
        Bbar_known = dfOIS_cf(idx_known) .* exp(-Z_known .* tau_cf(idx_known));
        knownPV    = sum(c_cf(idx_known) .* Bbar_known);
    end

    c_unk     = c_cf(idx_unknown);
    tau_unk   = tau_cf(idx_unknown);
    dfOIS_unk = dfOIS_cf(idx_unknown);

    if k == 1
        Z_prev = 0;
        w      = ones(size(tau_unk));
    else
        Z_prev   = Z(k-1);
        tau_prev = yearfrac(t0, T(k-1), 3);
        tau_k    = yearfrac(t0, T(k),   3);
        w        = (tau_unk - tau_prev) / (tau_k - tau_prev);
    end

    Z_unk     = @(z) Z_prev + (z - Z_prev) .* w;
    Bbar_unk  = @(z) dfOIS_unk .* exp(-Z_unk(z) .* tau_unk);
    unknownPV = @(z) sum(c_unk .* Bbar_unk(z));

    f    = @(z) knownPV + unknownPV(z) - target;
    Z(k) = fzero(f, [-0.01, 0.50]);
end

ZS.T     = T;
ZS.Z     = Z;
dfOIS_T  = interpolateOIS(t0, OIS.T, OIS.DF, T);
ZS.DFbar = dfOIS_T .* exp(-Z .* yearfrac(t0, T, 3));

tau = yearfrac(t0, T, 3);
Zbp = Z * 1e4;

fprintf('\nZ-spread bootstrap -- %s\n', name);
fprintf('  %-14s  %8s  %10s  %12s\n', 'Maturity', 'tau (y)', 'Z (bp)', 'Bbar');
fprintf('  %s\n', repmat('-', 1, 50));
for k = 1:N
    fprintf('  %-14s  %8.4f  %10.2f  %12.7f\n', ...
        datestr(T(k), 'dd-mmm-yyyy'), tau(k), Zbp(k), ZS.DFbar(k));
end

end