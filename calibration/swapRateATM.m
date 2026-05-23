function [S_atm, BPV, N_float] = swapRateATM(t_alpha, n_years, OIS, EUR6M)
%% SWAPRATEATM  Multicurve ATM forward swap rate.
%
%  Implements S_{alpha,omega}(t0)=N/BPV consistently with Section 2.3 of
%  Baviera (2019), for an Euribor 6m swap with annual fixed leg and
%  semi-annual floating leg.

t_fix     = t_alpha + (1:n_years)';
delta_fix = ones(n_years, 1);
B_fix     = getDF(OIS, t_fix);
BPV       = sum(delta_fix .* B_fix);

dt_float = 0.5;
t_float  = t_alpha + (0:2*n_years)' * dt_float;
PD_float = getPD(EUR6M, t_float);
B_float  = getDF(OIS, t_float);

PD_fwd_ratio = PD_float(1:end-1) ./ PD_float(2:end);
N_float      = sum((PD_fwd_ratio - 1) .* B_float(2:end));
S_atm        = N_float / BPV;
end
