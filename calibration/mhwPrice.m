function [price, xstar, S_atm] = mhwPrice(a, sigma, gamma, t_alpha, n_years, OIS, EUR6M)
%% MHWPRICE  MHW ATM CS receiver swaption price – Proposition 3.2, Baviera (2019).

%% 1. Discount factor at expiry
B_alpha = getDF(OIS, t_alpha);
n       = n_years;

%% 2. zeta_alpha (eq. 3.3)
zeta_a = mhwZeta(a, sigma, t_alpha);

%% 3. Fixed-leg DFs  (delta=1 yr per coupon period)
t_fix  = t_alpha + (1:n)';
B_fix  = getDF(OIS, t_fix) / B_alpha;
dlt_fx = ones(n, 1);

%% 4. Floating-leg dates and forward factors
t_fl   = t_alpha + (0:2*n)' * 0.5;
B_fl   = getDF(OIS, t_fl) / B_alpha;
PD_al  = getPD(EUR6M, t_alpha);
PD_fl  = getPD(EUR6M, t_fl) / PD_al;

% beta_i(t0)*B_{alpha,i+1}(t0) / PD_{alpha,i+1}(t0)  for i=0,...,2n-1
beta_B = PD_fl(1:end-1) .* B_fl(2:end) ./ PD_fl(2:end);

%% 5. Vol coefficients  (Lemma 3.1, eq. 3.5-3.7)
v_fix  = mhwV(a, zeta_a, t_fix - t_alpha);
v_fl   = mhwV(a, zeta_a, t_fl  - t_alpha);
xi_fix = (1 - gamma) * v_fix;
xi_fl  = (1 - gamma) * v_fl;
nu_fl  = v_fl(1:end-1) - gamma * v_fl(2:end);

%% 6. ATM swap rate S_atm  (eq. 2.14 / 2.15)
BPV_0  = sum(dlt_fx .* B_fix);
N_0    = beta_B(1) - B_fl(end) + sum(beta_B(2:end) - B_fl(2:end-1));
S_atm  = N_0 / BPV_0;

%% 7. S(x) function handle
S_fun = @(x) computeSwapRate(x, xi_fix, dlt_fx, B_fix, xi_fl, nu_fl, B_fl, beta_B);

%% 8. Find x*  (Lemma 3.2: unique zero of S(x) - S_atm)
%
%  Grid search over [-6, 6] to bracket the sign change, then fzero.
%  The fzero call is wrapped in its own try/catch so that a failure
%  (e.g. S_fun returns NaN at x=0) does not propagate upward as an
%  exception; instead we return price=NaN and let the caller handle it.

f_zero   = @(x) S_atm - S_fun(x);
x_grid   = linspace(-6, 6, 301);
f_grid   = arrayfun(f_zero, x_grid);
sign_chg = find(diff(sign(f_grid)) ~= 0, 1);

try
    if ~isempty(sign_chg)
        xstar = fzero(f_zero, [x_grid(sign_chg), x_grid(sign_chg+1)]);
    else
        % No sign change on grid: try searching from x=0
        xstar = fzero(f_zero, 0);
    end
catch
    % fzero failed: return NaN so swaptionObjective adds penalty
    price = NaN;
    return;
end

%% 9. Integration via composite Simpson on 400 subintervals
%
%  Each integrand evaluation uses integrandSafe which ALWAYS returns a
%  finite value (returns 0 for any degenerate input).
%
%  Using a fixed grid avoids the silent NaN propagation that can occur
%  inside MATLAB's adaptive integral().

x_lo = -8;
if xstar <= x_lo
    price = 0;
    return;
end

N_pts = 400;                                  % must be even
xv    = linspace(x_lo, xstar, N_pts+1);
dx    = (xstar - x_lo) / N_pts;

fv    = zeros(1, N_pts+1);
for k = 1:N_pts+1
    fv(k) = integrandSafe(xv(k), S_fun, S_atm, n_years);
end

% Composite Simpson weights: 1, 4, 2, 4, 2, ..., 4, 1
w            = 2*ones(1, N_pts+1);
w(2:2:end-1) = 4;
w(1)         = 1;
w(end)       = 1;

price_int = (dx/3) * (w * fv');
price     = B_alpha * price_int;
end


%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

function S = computeSwapRate(x, xi_fix, dlt_fx, B_fix, xi_fl, nu_fl, B_fl, beta_B)
% S(x) = N(x) / BPV(x),  formula (3.12) of Baviera (2019).
BPV      = sum(dlt_fx .* B_fix .* exp(-xi_fix*x - 0.5*xi_fix.^2));
term_b0  =  beta_B(1)  * exp(-nu_fl(1)*x    - 0.5*nu_fl(1)^2);
term_Bom = -B_fl(end)  * exp(-xi_fl(end)*x  - 0.5*xi_fl(end)^2);
if length(nu_fl) > 1
    term_spr = sum(beta_B(2:end)  .* exp(-nu_fl(2:end)*x   - 0.5*nu_fl(2:end).^2));
    term_dis = sum(B_fl(2:end-1)  .* exp(-xi_fl(2:end-1)*x - 0.5*xi_fl(2:end-1).^2));
else
    term_spr = 0;  term_dis = 0;
end
N = term_b0 + term_Bom + term_spr - term_dis;
if abs(BPV) < 1e-15,  S = NaN;  else,  S = N/BPV;  end
end


function val = integrandSafe(x, S_fun, S_atm, n_years)
% phi(x)*C(S(x))*[S_atm-S(x)]+  – returns 0 for any degenerate input.
try
    S_x = S_fun(x);
    if ~isfinite(S_x),       val=0; return; end
    payoff = S_atm - S_x;
    if payoff < 1e-15,       val=0; return; end
    if S_x <= -1 + 1e-10,   val=0; return; end
    C_x = cashAnnuity(S_x, n_years, 1);
    if ~isfinite(C_x)||C_x<=0, val=0; return; end
    val = normpdf(x) * C_x * payoff;
catch
    val = 0;
end
end
