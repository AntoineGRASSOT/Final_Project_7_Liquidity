function [S_atm, BPV, N_float] = swapRateATM(t_alpha, n_years, OIS, EUR6M)
%% SWAPRATEATM  Multicurve ATM forward swap rate.
%
%   S_{alpha,omega}(t0) = N_float / BPV  
%  for a Euribor 6m swap with EUR market conventions:
%    - fixed leg  : annual,       30/360   (basis=6)
%    - floating leg: semi-annual, ACT/360  (basis=2)
%
%  Payment dates are generated via addMonths (modified-following convention).
%
%  Inputs:
%    t_alpha  - swap start as year fraction ACT/365 from OIS.t0
%    n_years  - tenor in years (integer)
%    OIS      - struct: t0, T, DF   (bootstrap format)
%    EUR6M    - struct: t0, T, PD
%
%  Outputs:
%    S_atm   - ATM forward swap rate
%    BPV     - Basis Point Value (annuity), 30/360 day-count
%    N_float - floating leg numerator

t0            = OIS.t0;
expiry_y      = t_alpha;
expiry_months = round(12 * expiry_y);
d_alpha       = addMonths(t0, expiry_months);
T_alpha       = yearfrac(t0, d_alpha, 3);   % ACT/365

%% ── Fixed leg  ────────────────────────────────
d_fix     = addMonths(d_alpha, (1:n_years)' * 12);
t_fix     = yearfrac(t0, d_fix, 3);                   % ACT/365 for curve query

d_fix_all = [d_alpha; d_fix];
delta_fix = yearfrac(d_fix_all(1:end-1), d_fix_all(2:end), 6);  

B_fix = getDF(OIS, t_fix);
BPV   = sum(delta_fix .* B_fix);

%% ── Floating leg  ────────
d_float = addMonths(d_alpha, (0:2*n_years)' * 6);
t_float = yearfrac(t0, d_float, 3);

PD_float     = getPD(EUR6M, t_float);
B_float      = getDF(OIS,   t_float);

PD_fwd_ratio = PD_float(1:end-1) ./ PD_float(2:end);
N_float      = sum((PD_fwd_ratio - 1) .* B_float(2:end));

S_atm = N_float / BPV;
end
