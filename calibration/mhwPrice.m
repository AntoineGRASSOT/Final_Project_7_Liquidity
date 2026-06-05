function [price, xstar, S_atm] = mhwPrice(a, sigma, gamma, t_alpha, n_years, OIS, EUR6M)
%% MHWPRICE  MHW ATM CS receiver swaption price.
%
%  Euribor-6m CS receiver swaption with annual fixed leg and semi-annual floating leg.
%
%  Fixed leg:    annual,       30/360   (EUR market convention)
%  Floating leg: semi-annual,  ACT/360  (Euribor 6m convention)
%
%  Inputs:
%    a, sigma, gamma  - MHW model parameters 
%    t_alpha          - swaption expiry as year fraction ACT/365 from OIS.t0
%    n_years          - swap tenor in years 
%    OIS              - struct: t0, T, DF  
%    EUR6M            - struct: t0, T, PD
%
%  Outputs:
%    price   - CS receiver swaption price at t0
%    xstar   - critical value x* 
%    S_atm   - ATM forward swap rate


%% 1. ATM swap rate and schedule/curves
[S_atm, ~, ~, s] = swapRateATM(t_alpha, n_years, OIS, EUR6M);

T_alpha = s.T_alpha;
B_alpha = s.B_float(1);      
PD_al   = s.PD_float(1);    

zeta_a  = mhwZeta(a, sigma, T_alpha);

%% 2. Fixed leg (forward OIS DFs)
t_fix  = s.t_fix;
dlt_fx = s.delta_fix;                 % 30/360
B_fix  = s.B_fix / B_alpha;           % forward OIS DFs

%% 3. Floating leg (forward OIS / pseudo DFs)
t_fl   = s.t_float;
B_fl   = s.B_float  / B_alpha;        % forward OIS DFs
PD_fl  = s.PD_float / PD_al;          % forward pseudo-DFs

%% beta_B(k) = B_{alpha'k+1}(t0) * beta_k(t0)  
beta_B = (PD_fl(1:end-1) ./ PD_fl(2:end)) .* B_fl(2:end);   % 2n x 1

%% 4. Vol coefficients
%   v_{alpha',i} = zeta_a * (1-exp(-a*tau_i))/a
%   xi_{alpha',i} = (1-gamma) * v_{alpha',i}  
%   nu_{alpha',i} = v_{alpha',i} - gamma * v_{alpha',i+1}  

tau_fix = t_fix - T_alpha;            % year fracs from expiry  (n x 1)
tau_fl  = t_fl  - T_alpha;            % year fracs from expiry  (2n+1 x 1)

v_fix  = mhwV(a, zeta_a, tau_fix);
v_fl   = mhwV(a, zeta_a, tau_fl);

xi_fix = (1 - gamma) * v_fix;
xi_fl  = (1 - gamma) * v_fl;
nu_fl  = v_fl(1:end-1) - gamma * v_fl(2:end);


%% 5. S(x) function 
S_fun = @(x) computeSwapRate(x, xi_fix, dlt_fx, B_fix, xi_fl, nu_fl, B_fl, beta_B);

%% 6. Find x*  (unique zero of S_atm - S(x))
f_zero   = @(x) S_atm - S_fun(x);
x_grid   = linspace(-6, 6, 301);
f_grid   = f_zero(x_grid);                           
sign_chg = find(diff(sign(f_grid)) ~= 0, 1);

try
    if ~isempty(sign_chg)
        xstar = fzero(f_zero, [x_grid(sign_chg), x_grid(sign_chg + 1)]);
    else
        xstar = fzero(f_zero, 0);
    end
catch
    price = NaN;
    return;
end

%% 7. Integration - composite Simpson on 400 subintervals
x_lo = -8;
if xstar <= x_lo
    price = 0;
    return;
end

N_pts = 400;                                           
xv    = linspace(x_lo, xstar, N_pts + 1);             
dx    = (xstar - x_lo) / N_pts;

fv = integrandVec(xv, S_fun, S_atm, n_years);         

% Composite Simpson weights: 1, 4, 2, 4, 2, ..., 4, 1
w            = 2 * ones(1, N_pts + 1);
w(2:2:end-1) = 4;
w(1)         = 1;
w(end)       = 1;

price_int = (dx / 3) * (w * fv');
price     = B_alpha * price_int;
end


%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

function S = computeSwapRate(x, xi_fix, dlt_fx, B_fix, xi_fl, nu_fl, B_fl, beta_B)
%% S(x) = N(x) / BPV(x)  
%
%  BPV(x) = sum_j dlt_j * B_fix_j * exp(-xi_fix_j * x - xi^2/2)
%  N(x)   = sum_k beta_B_k        * exp(-nu_k      * x - nu^2/2)
%           -sum_k B_fl(k+1)      * exp(-xi_fl(k+1)* x - xi^2/2)

x = x(:)';   % ensure 1xM row for broadcasting

BPV = sum(dlt_fx .* B_fix .* exp(-xi_fix .* x - 0.5 * xi_fix.^2), 1);

% First sum: beta_B uses end-of-period B_fl(k+1), paired with nu_fl(k)
N1  = sum(beta_B .* exp(-nu_fl .* x - 0.5 * nu_fl.^2), 1);

% Second sum: B_fl(2:end) = B_{alpha',k+1} for k=1..2n, paired with xi_fl(2:end)
N2  = sum(B_fl(2:end) .* exp(-xi_fl(2:end) .* x - 0.5 * xi_fl(2:end).^2), 1);

N = N1 - N2;
S = N ./ BPV;
S(abs(BPV) < 1e-15) = NaN;
end


function fv = integrandVec(xv, S_fun, S_atm, n_years)
%% phi(x) * C(S(x)) * [S_atm - S(x)]+  

xv  = xv(:)';
S_x = S_fun(xv);
fv  = zeros(1, numel(xv));

valid = isfinite(S_x) & (S_atm - S_x > 1e-15) & (S_x > -1 + 1e-10);
if ~any(valid), return; end

C_x   = cashAnnuity(S_x(valid), n_years, 1);   
fin_C = isfinite(C_x) & (C_x > 0);

mask        = valid;
mask(valid) = fin_C;

fv(mask) = normpdf(xv(mask)) .* C_x(fin_C) .* (S_atm - S_x(mask));
end
