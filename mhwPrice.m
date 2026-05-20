function [price, xstar, S_atm] = mhwPrice(a, sigma, gamma, t_alpha, n_years, OIS, EUR6M)
%% MHWPRICE  MHW price of an ATM CS receiver swaption (Proposition 3.2)
%
%  Implements the closed formula (3.13) in [B2019]:
%
%    R^{C,MHW}_{αω}(t0) = B(t0,tα) * ∫_{-∞}^{x*} φ(x) * C_{αω}(S(x)) * [K-S(x)]⁺ dx
%
%  where phi(x) = exp(-x^2/2)/sqrt(2*pi) is the standard normal density and x is the
%  standard Gaussian random variable defined in (3.2).
%
%  --------------------------------------------------------------------------
%  MHW MODEL INGREDIENTS (Lemma 3.1, eq. 3.3-3.7):
%
%  (1) zeta_alpha^2 = sigma^2(1 - exp(-2a*t_alpha)) / (2a) [factor variance]
%      zeta_alpha = sqrt(zeta_alpha^2)
%
%  (2) For each floating date t_i or fixed date t_j:
%      v_{alpha,i} = zeta_alpha * (1 - exp(-a*(t_i-t_alpha))) / a [integrated HJM coefficient]
%
%  (3) Volatility coefficients (eq. 3.6-3.7):
%      xi_{alpha,i} = (1-gamma) * v_{alpha,i} [for forward discounts]
%      nu_{alpha,i} = v_{alpha,i} - gamma * v_{alpha,i+1} [for forward spreads, extended vol]
%
%  (4) Swap rate S(x) from formula (3.12):
%      S(x) = N(x) / BPV(x)
%
%      BPV(x) = Σ_{j=1}^n δⱼ B_{αj}(t0) exp(-ς_j^fix * x - ½ς_j^fix²)
%
%      N(x) = β_0(t0) * exp(-ν_0*x - ½ν_0²)                           [ι=0]
%            - B_{α,2n}^float(t0) * exp(-ς_{2n}*x - ½ς_{2n}²)         [-B_alpha,omega term]
%            + Σ_{ι=1}^{2n-1} β_ι(t0)B_{αι}^float * exp(-ν_ι*x-½ν_ι²) [spread]
%            - Σ_{ι=1}^{2n-1} B_{αι}^float(t0) * exp(-ς_ι*x - ½ς_ι²) [discount]
%
%      where:
%        B_{alpha,j}(t0) = B(t0,t_j)/B(t0,t_alpha) [fixed-leg forward discount]
%        B_{alpha,i}(t0) = B(t0,t'_i)/B(t0,t_alpha) [floating-date forward discount]
%        β_ι(t0)*B_{αι}(t0) = PD(t0,t'ι)*B(t0,t'ι+1) / [PD(t0,t'ι+1)*B(t0,tα)]
%
%  (5) x* is the unique zero of f(x) = S(x) - S_ATM (Lemma 3.2):
%      There is a unique x* such that S(x*) = S_ATM = K (ATM swaption).
%      NOTE: f is not monotone in general (Fig. 2 of the paper),
%      but uniqueness is guaranteed by Lemma 3.2.
%
%  (6) Numerical integral (Gauss-Hermite or fine grid):
%      R^{C,MHW} = B(t0,tα) * ∫_{-∞}^{x*} φ(x) C(S(x)) (S_ATM - S(x)) dx
%
%  --------------------------------------------------------------------------
%  NOTE ON FLOATING LEG vs FIXED LEG:
%  - Fixed leg: annual payments at t_alpha+1, t_alpha+2, ..., t_alpha+n
%    delta_fix approximately 1 year (30/360)
%  - Floating leg: semi-annual resets/payments (Euribor 6m)
%    Dates: t'_i = t_alpha + i*0.5 for i = 0,...,2n
%    i=0: t'_0 = t_alpha (first reset date = swaption expiry)
%    i=2n: t'_{2n} = t_alpha + n = t_omega (last payment date = swap maturity)
%
%  --------------------------------------------------------------------------
%  Input:
%    a       - mean reversion speed, e.g. 0.1294
%    sigma   - pseudo-discount volatility scale, e.g. 0.0126
%    gamma   - fraction of volatility assigned to the spread [0,1]
%    t_alpha - swaption expiry in years from settlement, e.g. 1.0, 2.0, ...
%    n_years - swap tenor in years, e.g. 9, 8, ..., 1
%    OIS     - OIS curve structure (fields .T, .DF)
%    EUR6M   - Euribor 6m pseudo-discount curve structure (fields .T, .PD)
%
%  Output:
%    price  - ATM CS receiver swaption price per unit notional
%    xstar  - critical value x* (Lemma 3.2)
%    S_atm  - current ATM swap rate S_{alpha,omega}(t0)
%  --------------------------------------------------------------------------

%% -- Step 1: Basic quantities ---------------------------------------------
B_alpha = getDF(OIS, t_alpha);       % B(t0, tα)
n       = n_years;                   % tenor in years
dt_fl   = 0.5;                       % floating-leg step (6m)

%% -- Step 2: ζ_α (eq. 3.3) -----------------------------------------------
zeta_a = mhwZeta(a, sigma, t_alpha);

%% -- Step 3: Dates and forward quantities ------------------------------------

% Fixed leg: t_j = t_alpha + j years, j = 1,...,n
t_fix  = t_alpha + (1:n)';
B_fix  = getDF(OIS, t_fix) / B_alpha;     % B_{αj}(t0) = B(t0,tj)/B(t0,tα)
dlt_fx = ones(n,1);                        % delta_fix approximately 1 year (30/360)

% Floating leg: t'_i = t_alpha + i*0.5, i = 0,...,2n
t_fl   = t_alpha + (0:2*n)' * dt_fl;      % 2n+1 dates
B_fl   = getDF(OIS, t_fl) / B_alpha;      % B_{αι}^float(t0)
PD_fl  = getPD(EUR6M, t_fl) / getPD(EUR6M, t_alpha);  % PD_{αι}(t0) = PD(t0,t'ι)/PD(t0,tα)

% β_ι(t0) * B_{αι}(t0) for each floating period [ι=0,...,2n-1]:
%   = PD(t0,t'ι) * B(t0,t'ι+1) / [PD(t0,t'ι+1) * B(t0,tα)]
%   = PD_{αι} * (B_fl(ι+1) * B_alpha) / (PD_fl(ι+1) * B_alpha)
%   = PD_{αι}(t0) * B_{α,ι+1}(t0) / PD_{α,ι+1}(t0)
beta_B = PD_fl(1:end-1) .* B_fl(2:end) ./ PD_fl(2:end);   % 2n values [ι=0..2n-1]

%% -- Step 4: Coefficients v, xi, nu (eq. 3.6-3.7) -------------------------

% v_{αι} = ζ_α * (1 - e^{-a*(tι-tα)}) / a
v_fix  = mhwV(a, zeta_a, t_fix  - t_alpha);  % for fixed-leg dates
v_fl   = mhwV(a, zeta_a, t_fl   - t_alpha);  % for floating-leg dates

% ς_{αj}^fix = (1-γ) * v_j^fix
xi_fix = (1 - gamma) * v_fix;    % n values [j=1..n]

% ς_{αι}^float = (1-γ) * v_ι^float  [ι=0..2n]
xi_fl  = (1 - gamma) * v_fl;     % 2n+1 values

% ν_{αι} = v_{αι} - γ * v_{α,ι+1}  [ι=0..2n-1]
nu_fl  = v_fl(1:end-1) - gamma * v_fl(2:end);  % 2n values

%% -- Step 5: ATM swap rate S_atm ----------------------------------------
% S_atm = N(t0) / BPV(t0)
% Using forward quantities without the exponential Ito adjustment:

BPV_0 = sum(dlt_fx .* B_fix);

% N_0 through relation (2.14):
% N = 1 - B_{α,ω} + Σ_ι B_{αι}(β_ι-1)
% Equivalently:
% N = -B_{α,2n} + β_0*B_{α,0} + Σ_{ι=1}^{2n-1}(β_ι B_{αι} - B_{αι})
N_0 = beta_B(1) ...                            % β_0 * B_{α,0} = β_0 * 1
    - B_fl(end) ...                             % -B_{α,2n}
    + sum(beta_B(2:end) - B_fl(2:end-1));       % Σ_{ι=1}^{2n-1}(β_ι B_ι - B_ι)

S_atm = N_0 / BPV_0;

%% -- Step 6: Function S(x) -----------------------------------------------
%  Implements formula (3.12) in vectorized form over x

S_fun = @(x) computeSwapRate(x, ...
    xi_fix, dlt_fx, B_fix, ...      % fixed leg
    xi_fl, nu_fl, B_fl, beta_B);   % floating leg

%% -- Step 7: Find x* (Lemma 3.2) ----------------------------------------
%  x* is the value such that S(x*) = S_atm (ATM condition)
%  f(x) = S_atm - S(x) has a unique zero (Lemma 3.2)
%
%  NOTE: f(x) is NOT monotone in general (see Fig. 2 of the paper).
%  However, Lemma 3.2 guarantees uniqueness. We use fzero with initial point 0.

f_zero = @(x) S_atm - S_fun(x);

% Find an interval containing the zero
x_grid = linspace(-5, 5, 201);
f_grid = arrayfun(f_zero, x_grid);
% Look for the first sign change
sign_chg = find(diff(sign(f_grid)) ~= 0, 1);
if ~isempty(sign_chg)
    x_lo = x_grid(sign_chg);
    x_hi = x_grid(sign_chg + 1);
    try
        xstar = fzero(f_zero, [x_lo, x_hi]);
    catch
        xstar = fzero(f_zero, 0);
    end
else
    % Fallback: initial point 0
    xstar = fzero(f_zero, 0);
end

%% -- Step 8: Numerical integral (Proposition 3.2) -----------------------
%  R^{C,MHW} = B(t0,tα) * ∫_{-∞}^{x*} φ(x) * C(S(x)) * [S_atm - S(x)] dx
%
%  Uses Gaussian integration on [-n_sigma, x*], since x is already N(0,1).
%  Since S(x) is smooth, quadrature is efficient.

n_sigma = 8;   % integrate from -8; contribution below -8 is negligible
x_lo_int = -n_sigma;

if xstar <= x_lo_int
    % very negative x*: price is almost zero (deep OTM)
    price = 0;
    return;
end

% Integration with MATLAB integral (adaptive, high precision)
integrand = @(x) arrayfun(@(xi) integrandFun(xi, S_fun, S_atm, n_years), x);
price_int = integral(integrand, x_lo_int, xstar, ...
                     'RelTol', 1e-6, 'AbsTol', 1e-10);
price = B_alpha * price_int;
end


%% =========================================================================
%%  LOCAL HELPER FUNCTIONS
%% =========================================================================

function S = computeSwapRate(x, xi_fix, dlt_fx, B_fix, xi_fl, nu_fl, B_fl, beta_B)
%  Computes S(x) = N(x)/BPV(x) for a scalar x.
%  Implements formula (3.12) in [B2019].
%
%  BPV(x) = Σ_{j=1}^n δⱼ B_{αj}(t0) exp(-ς_j^fix * x - ½ς_j^fix²)
%
%  N(x) = β_0 B_{α,0}(t0) exp(-ν_0*x - ½ν_0²)          [ι=0: B_{α,0}=1]
%        - B_{α,2n}(t0) exp(-ς_{2n}*x - ½ς_{2n}²)       [-B_{alpha,omega} term]
%        + Σ_{ι=1}^{2n-1} βB_ι exp(-ν_ι*x - ½ν_ι²)      [spread ι=1..2n-1]
%        - Σ_{ι=1}^{2n-1} B_ι(t0) exp(-ς_ι*x - ½ς_ι²)  [discount ι=1..2n-1]

n2 = length(nu_fl);   % = 2n

% BPV(x): sum over the fixed leg j=1..n
BPV = sum(dlt_fx .* B_fix .* exp(-xi_fix*x - 0.5*xi_fix.^2));

% Term β_0 * B_{α,0} (ι=0: B_{α,0}=1, ν_0=nu_fl(1))
term_beta0 = beta_B(1) * exp(-nu_fl(1)*x - 0.5*nu_fl(1)^2);

% Term -B_{alpha,2n}; this is B_fl(end)
% ς_{2n} = xi_fl(end)
term_Bomega = -B_fl(end) * exp(-xi_fl(end)*x - 0.5*xi_fl(end)^2);

% Sums for i=1..2n-1
if n2 > 1
    % MATLAB indices: ι=1..2n-1 => indices 2..2n in vectors (1-indexed)
    term_spread  = sum(beta_B(2:end)   .* exp(-nu_fl(2:end)*x   - 0.5*nu_fl(2:end).^2));
    term_disc    = sum(B_fl(2:end-1)   .* exp(-xi_fl(2:end-1)*x  - 0.5*xi_fl(2:end-1).^2));
else
    term_spread = 0;
    term_disc   = 0;
end

N = term_beta0 + term_Bomega + term_spread - term_disc;

if abs(BPV) < 1e-15
    S = NaN;
else
    S = N / BPV;
end
end


function val = integrandFun(x, S_fun, S_atm, n_years)
%  Integrand for Proposition 3.2:
%  val = φ(x) * C(S(x)) * [S_atm - S(x)]⁺
%
%  Where C(S) = cashAnnuity(S, n_years, 1)

S_x = S_fun(x);
payoff = max(S_atm - S_x, 0);

if payoff < 1e-15
    val = 0;
    return;
end

C_x = cashAnnuity(S_x, n_years, 1);   % C_{αω}(S(x)), m=1 EUR
phi = normpdf(x);                       % φ(x) = exp(-x²/2)/√(2π)

val = phi * C_x * payoff;
end
