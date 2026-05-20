function [MHWparams, calibRes] = calibrateMHW(swaptVol, OIS, EUR6M, params)
%% CALIBRATEMHW  MHW calibration on 9 diagonal co-terminal 10y ATM CS swaptions
%
%  Implements the calibration procedure described in Section 4 of [B2019]:
%
%  "Second, we calibrate the three MHW parameters p := (a, σ, γ) with
%   European CS swaptions versus Euribor 6m on the 10y-diagonal
%   (i.e. considering M=9 swaptions 1y9y, 2y8y,...,9y1y)"
%
%  OBJECTIVE FUNCTION (eq. 4.1 in [B2019]):
%
%    Err²(p) = Σ_{i=1}^M [R^{C,MHW}_i(p;t0) - R^{C,mkt}_i(t0)]²
%
%  where:
%    R^{C,MHW}_i = MHW price from mhwPrice.m (Proposition 3.2)
%    R^{C,mkt}_i = market price from bachelierCS (formula B.4)
%                 = B(t0,t_alpha_i) * C(Si) * sigma^mkt_i * sqrt(t_alpha_i/(2*pi))   [ATM]
%
%  Prices, not volatilities, are compared to remain consistent with the
%  metric used in the paper. MHW implied volatilities are derived only for
%  post-calibration analysis.
%
%  PARAMETERS TO CALIBRATE:
%    a     ∈ (0, 1]      mean reversion speed
%    sigma in (0, 0.10]  pseudo-discount volatility scale
%    gamma in [0, 1]     weight of the spread dynamics
%
%  EXPECTED VALUES (Table 4, [B2019]):
%    a = 12.94%,  σ = 1.26%,  γ = 0.07%
%
%  --------------------------------------------------------------------------
%  NOTE ON gamma (Section 4 of [B2019]):
%  "The dependence of Err on γ is less pronounced compared to a and σ;
%   even if minimum values for Err are achieved for very low values of γ,
%   differences in squared distance are very small when increasing γ."
%  This explains why gamma is close to 0 in the S0 limiting case (deterministic spread).
%
%  --------------------------------------------------------------------------
%  Input:
%    swaptVol - swaption volatility structure with fields:
%               .expiry  - expiry in years [1,2,...,9] (M values)
%               .tenor   - tenor in years [9,8,...,1] (M values)
%               .vol     - normal volatility in DECIMAL form (M values)
%               .volBps  - normal volatility in bps (M values)
%    OIS      - OIS curve structure (.T, .DF)
%    EUR6M    - pseudo-discount curve structure (.T, .PD)
%    params   - parameter structure, kept for extensibility
%
%  Output:
%    MHWparams - structure with fields:
%                .a, .sigma, .gamma     (optimal parameters)
%                .a_pct, .sigma_pct, .gamma_pct (percent values)
%    calibRes  - structure with fields:
%                .price_mkt    - market prices R^{C,mkt}_i
%                .price_model  - MHW prices R^{C,MHW}_i at the optimum
%                .vol_mkt_bps  - market volatilities in bps
%                .vol_impl_bps - MHW implied volatilities in bps
%                .error_bps    - volatility difference (model-market) in bps
%                .RMSE_price   - price RMSE (objective function)
%                .RMSE_vol_bps - volatility RMSE in bps
%                .expiry, .tenor
%                .exitflag, .output (from fmincon)
%  --------------------------------------------------------------------------

fprintf('>>> MHW calibration on %d diagonal ATM swaptions...\n', ...
        length(swaptVol.expiry));

M = length(swaptVol.expiry);
expiry_y = swaptVol.expiry(:);   % [1,2,...,9] years
tenor_y  = swaptVol.tenor(:);   % [9,8,...,1] years
vol_mkt  = swaptVol.vol(:);      % decimal form

%% -- Step 1: Market prices R^{C,mkt}_i --------------------------------------
%  ATM Bachelier formula (Appendix B, eq. B.4 simplified for ATM):
%    R^{C,mkt}_{ATM} = B(t0,tα) * C(S0) * σ_mkt * sqrt(tα/(2π))
%
%  NOTE: the square-root term includes the 1/(2*pi) factor from the ATM formula.

price_mkt = zeros(M, 1);
S_atm_vec = zeros(M, 1);
B_alpha_v = zeros(M, 1);
C_S0_vec  = zeros(M, 1);

for i = 1:M
    t_a  = expiry_y(i);
    n_y  = tenor_y(i);
    sig  = vol_mkt(i);

    B_a  = getDF(OIS, t_a);
    [S0, ~, ~] = swapRateATM(t_a, n_y, OIS, EUR6M);
    C0 = cashAnnuity(S0, n_y, 1);

    % ATM Bachelier CS receiver swaption price (eq. B.3/B.4 for ATM)
    % Calls bachelierCS.m (K=S0 for ATM)
    price_mkt(i) = bachelierCS(S0, S0, t_a, B_a, C0, sig);

    B_alpha_v(i) = B_a;
    S_atm_vec(i) = S0;
    C_S0_vec(i)  = C0;
end

fprintf('   Market prices computed (range: [%.6f, %.6f])\n', ...
        min(price_mkt), max(price_mkt));

%% -- Step 2: Objective function ------------------------------------------
%  Err²(p) = Σ [R^{C,MHW}_i(p) - R^{C,mkt}_i]²

objFun = @(p) swaptionErr2(p(1), p(2), p(3), ...
                           expiry_y, tenor_y, OIS, EUR6M, price_mkt);

%% -- Step 3: Optimization with fmincon -----------------------------------
%  Physical bounds:
%    a     ∈ (0, 1]
%    σ     ∈ (0, 0.10]
%    γ     ∈ [0, 1]     (γ=0 = S0 hypothesis, γ=1 = S1 hypothesis)
%
%  Multiple starting points for robustness. The paper states that the solution
%  is stable for a broad class of initial points (Section 4, [B2019]).

p0_list = [
    0.10,  0.01,   0.001;   % close to the expected solution
    0.05,  0.015,  0.0005;
    0.20,  0.008,  0.002;
    0.15,  0.012,  0.0;     % γ=0 (S0 hypothesis)
];

lb = [1e-4,  1e-4,  0.0 ];
ub = [1.0,   0.10,  1.0 ];

opt = optimoptions('fmincon', ...
    'Algorithm',             'interior-point', ...
    'Display',               'off', ...
    'MaxIterations',         1000, ...
    'MaxFunctionEvaluations',5000, ...
    'FunctionTolerance',     1e-12, ...
    'StepTolerance',         1e-12, ...
    'OptimalityTolerance',   1e-12);

best_fval = Inf;
best_p    = p0_list(1,:);
best_ef   = -1;
best_out  = [];

for k = 1:size(p0_list, 1)
    p0_k = p0_list(k,:);
    try
        [p_k, fval_k, ef_k, out_k] = fmincon(objFun, p0_k, [], [], [], [], lb, ub, [], opt);
        fprintf('   Start %d: a=%.4f, σ=%.4f, γ=%.5f => Err=%.4e (exitflag=%d)\n', ...
                k, p_k(1), p_k(2), p_k(3), sqrt(fval_k), ef_k);
        if fval_k < best_fval
            best_fval = fval_k;
            best_p    = p_k;
            best_ef   = ef_k;
            best_out  = out_k;
        end
    catch ME
        fprintf('   Start %d: error (%s)\n', k, ME.message);
    end
end

a_opt     = best_p(1);
sigma_opt = best_p(2);
gamma_opt = best_p(3);

%% -- Step 4: MHW prices and volatilities with optimal parameters --------------------
price_model  = zeros(M, 1);
vol_impl_bps = zeros(M, 1);

for i = 1:M
    t_a = expiry_y(i);
    n_y = tenor_y(i);

    % MHW price
    price_model(i) = mhwPrice(a_opt, sigma_opt, gamma_opt, t_a, n_y, OIS, EUR6M);

    % MHW implied volatility (Bachelier inversion)
    price_i = price_model(i);
    B_a     = B_alpha_v(i);
    C0      = C_S0_vec(i);
    vol_impl_bps(i) = bachelierCS_implVol(S_atm_vec(i), S_atm_vec(i), t_a, B_a, C0, price_i) * 1e4;
end

vol_mkt_bps = vol_mkt * 1e4;
error_bps   = vol_impl_bps - vol_mkt_bps;
RMSE_price  = sqrt(best_fval / M);
RMSE_vol    = sqrt(mean(error_bps.^2));

%% -- Output ---------------------------------------------------------------
MHWparams.a         = a_opt;
MHWparams.sigma     = sigma_opt;
MHWparams.gamma     = gamma_opt;
MHWparams.a_pct     = a_opt     * 100;
MHWparams.sigma_pct = sigma_opt * 100;
MHWparams.gamma_pct = gamma_opt * 100;

calibRes.price_mkt    = price_mkt;
calibRes.price_model  = price_model;
calibRes.vol_mkt_bps  = vol_mkt_bps;
calibRes.vol_impl_bps = vol_impl_bps;
calibRes.error_bps    = error_bps;
calibRes.RMSE_price   = RMSE_price;
calibRes.RMSE_vol_bps = RMSE_vol;
calibRes.expiry       = expiry_y;
calibRes.tenor        = tenor_y;
calibRes.exitflag     = best_ef;
calibRes.output       = best_out;

%% -- Print results -------------------------------------------------------
fprintf('\n--- MHW CALIBRATION RESULTS ---\n');
fprintf('  a     = %8.4f%%  (expected: 12.94%%)\n', a_opt*100);
fprintf('  sigma = %8.4f%%  (expected:  1.26%%)\n', sigma_opt*100);
fprintf('  gamma = %8.5f%%  (expected:  0.07%%)\n', gamma_opt*100);
fprintf('  price RMSE    = %.4e\n', RMSE_price);
fprintf('  vol RMSE      = %.3f bps\n', RMSE_vol);
fprintf('\n  %-8s %-8s %-12s %-12s %-10s\n', ...
        'Expiry', 'Tenor', 'Mkt vol(bps)', 'Vol MHW(bps)', 'Err(bps)');
fprintf('  %s\n', repmat('-', 1, 55));
for i = 1:M
    fprintf('  %dy       %dy       %10.2f   %10.2f   %+8.3f\n', ...
            expiry_y(i), tenor_y(i), vol_mkt_bps(i), vol_impl_bps(i), error_bps(i));
end
fprintf('\n');
end


%% =========================================================================
%%  LOCAL HELPER FUNCTIONS
%% =========================================================================

function err2 = swaptionErr2(a, sigma, gamma, expiry_y, tenor_y, OIS, EUR6M, price_mkt)
%  Objective function: sum of squared PRICE differences (eq. 4.1)
M    = length(expiry_y);
err2 = 0;
for i = 1:M
    try
        p_i   = mhwPrice(a, sigma, gamma, expiry_y(i), tenor_y(i), OIS, EUR6M);
        err2  = err2 + (p_i - price_mkt(i))^2;
    catch
        err2 = err2 + 1e6;  % penalize failed evaluations
    end
end
end
