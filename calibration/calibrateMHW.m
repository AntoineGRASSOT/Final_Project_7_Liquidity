function [a_hat, sigma_hat, gamma_hat, calibRes] = calibrateMHW(data, dates, dates_OIS, df_OIS, E6m_dates, E6m_df)
%% CALIBRATEMHW  MHW calibration – 9 diagonal ATM CS swaptions.
%
%  Minimises:  Err^2 = sum_{i=1}^9 [R^{C,MHW}_i - R^{C,mkt}_i]^2
%  over p = (a, sigma, gamma) on co-terminal 10y diagonal:
%  1y9y, 2y8y, ..., 9y1y.
%
%  Inputs:
%    data      - readExcelData_OIS output  
%    dates     - readExcelData_OIS output  
%    dates_OIS - OIS serial dates    from bootstrapOIS
%    df_OIS    - OIS discount factors from bootstrapOIS
%    E6m_dates - EUR6M serial dates   from bootstrapEuribor6m
%    E6m_df    - EUR6M pseudo-DFs     from bootstrapEuribor6m
%
%  Outputs:
%    a_hat, sigma_hat, gamma_hat  (ref: 12.94%, 1.26%, 0.07%)
%    calibRes  - diagnostics struct

%% ── 0. Curve structs ─────────────────────────────────────────────────────
VERBOSE = false; %pre-check flag section 4. Could be used for the debug

t0       = dates.settlement;
OIS.t0   = t0;          
OIS.T    = yearfrac(t0, dates_OIS(:), 3);   
OIS.DF   = df_OIS(:);

EUR6M.t0 = t0;
EUR6M.T  = yearfrac(t0, E6m_dates(:), 3);
EUR6M.PD = E6m_df(:); 

fprintf('\n>>> MHW calibration: 9 diagonal ATM swaptions\n');

%% ── 1. Swaption grid ─────────────────────────────────────────────────────
M        = 9;
expiry_y = (1:M)';      % 1y … 9y
tenor_y  = (M:-1:1)';   % 9y … 1y

%% ── 2. Market vols  ──────────────────────────────────────────────────────

vols_all = data.vols(:);
idx0     = find(isfinite(vols_all), 1, 'first');   

if isempty(idx0)
    error('calibrateMHW: data.vols has no finite values.');
end
if numel(vols_all) < idx0 + M - 1
    error('calibrateMHW: fewer than %d finite vols available in data.vols.', M);
end

vol_bps = vols_all(idx0 : idx0+M-1);   
vol_dec = vol_bps / 1e4;               

fprintf('   Vol index: data.vols(%d:%d)\n', idx0, idx0+M-1);
fprintf('   Vols (bps):'); fprintf(' %6.2f', vol_bps); fprintf('\n');

%% ── 3. Market prices  R^{C,mkt}  (Bachelier ATM) ───────────────
p_mkt    = zeros(M,1);
S_vec    = zeros(M,1);
B_vec    = zeros(M,1);
C_vec    = zeros(M,1);

for i = 1:M
    T_exp    = yearfrac(t0, addMonths(t0, 12 * expiry_y(i)), 3);
    B_a      = getDF(OIS, T_exp);
    [S0,~,~] = swapRateATM(expiry_y(i), tenor_y(i), OIS, EUR6M);
    C0       = cashAnnuity(S0, tenor_y(i), 1);
    p_mkt(i) = bachelierCS(S0, S0, T_exp, B_a, C0, vol_dec(i));
    B_vec(i) = B_a;  S_vec(i) = S0;  C_vec(i) = C0;
end
fprintf('   Market prices: [%.4e , %.4e]\n', min(p_mkt), max(p_mkt));

%% ── 4. Objective ─────────────────────────────────────────────────────────
objFun = @(p) swaptionObj(p, expiry_y, tenor_y, OIS, EUR6M, p_mkt);

% Verbose pre-check: show per-swaption results at p0(1)
if VERBOSE
    fprintf('\n   Per-swaption check at p0=[0.10, 0.01, 0.001]:\n');
    a0=0.10; s0=0.01; g0=0.001;
    for i = 1:M
        try
            pv = mhwPrice(a0,s0,g0,expiry_y(i),tenor_y(i),OIS,EUR6M);
            sq = (pv - p_mkt(i))^2;
            fprintf('   %dy%dy: price=%.4e  pmkt=%.4e  sq=%.2e\n', ...
                    expiry_y(i),tenor_y(i), pv, p_mkt(i), sq);
        catch ME
            fprintf('   %dy%dy: ERROR -> %s\n', expiry_y(i),tenor_y(i), ME.message);
        end
    end
    fprintf('   Objective at p0(1): %.6e\n\n', objFun([a0,s0,g0]));
end

%% ── 5. fmincon ───────────────────────────────────────────────────────────
lb  = [1e-4, 1e-4, 0.0];
ub  = [1.0,  0.10, 1.0];
p0s = [0.10, 0.010, 0.001;
       0.05, 0.015, 0.0005;
       0.20, 0.008, 0.002;
       0.15, 0.012, 0.0];

opt = optimoptions('fmincon','Algorithm','interior-point','Display','off', ...
    'MaxIterations',2000,'MaxFunctionEvaluations',10000, ...
    'FunctionTolerance',1e-14,'StepTolerance',1e-12,'OptimalityTolerance',1e-12);

best_val = Inf;  best_p = p0s(1,:);  best_ef = -1;
for k = 1:size(p0s,1)
    try
        [pk,fk,efk] = fmincon(objFun, p0s(k,:), [],[],[],[], lb, ub, [], opt);
        fprintf('   Start %d: a=%7.4f%%  sig=%6.4f%%  gam=%8.5f%%  Err=%.3e  (flag=%d)\n', ...
                k, pk(1)*100, pk(2)*100, pk(3)*100, sqrt(fk), efk);
        if fk < best_val,  best_val=fk;  best_p=pk;  best_ef=efk;  end
    catch ME
        fprintf('   Start %d: fmincon error -> %s\n', k, ME.message);
    end
end

a_hat     = best_p(1);
sigma_hat = best_p(2);
gamma_hat = best_p(3);

%% ── 6. Model prices and implied vols at optimum ──────────────────────────
p_mod   = zeros(M,1);
vi_bps  = zeros(M,1);          % model implied vols in bps
for i = 1:M
    T_exp    = yearfrac(t0, addMonths(t0, 12 * expiry_y(i)), 3);
    p_mod(i) = mhwPrice(a_hat, sigma_hat, gamma_hat, ...
                         expiry_y(i), tenor_y(i), OIS, EUR6M);
    vi_bps(i) = bachelierCS_implVol(S_vec(i), S_vec(i), T_exp, ...
                    B_vec(i), C_vec(i), p_mod(i)) * 1e4;
end

err_bps    = vi_bps - vol_bps;      % in bps
RMSE_price = sqrt(best_val / M);
RMSE_vol   = sqrt(mean(err_bps.^2));

%% ── 6b. Calibration plot ────────────────────────────────────────────────
labels = arrayfun(@(e,t) sprintf('%dy%dy',e,t), expiry_y, tenor_y, ...
                  'UniformOutput', false);

figure('Name','MHW calibration','Color','w','Position',[100 100 760 540]);

subplot(2,1,1);
plot(expiry_y, vol_bps, 'o-', 'LineWidth',1.4, 'MarkerSize',6); hold on;
plot(expiry_y, vi_bps,  's--','LineWidth',1.4, 'MarkerSize',6);
grid on; box on; xticks(expiry_y); xticklabels(labels);
ylabel('Normal vol (bps)');
legend('Market','MHW','Location','best');
title(sprintf('Co-terminal 10y diagonal — Vol RMSE = %.2f bps', RMSE_vol));

subplot(2,1,2);
bar(expiry_y, err_bps, 0.6);
grid on; box on; xticks(expiry_y); xticklabels(labels);
ylabel('Error (bps)'); xlabel('Swaption');
title('MHW - Market implied vol');

%% ── 7. Output struct ─────────────────────────────────────────────────────
calibRes.price_mkt    = p_mkt;
calibRes.price_model  = p_mod;
calibRes.vol_mkt_bps  = vol_bps;
calibRes.vol_impl_bps = vi_bps;
calibRes.error_bps    = err_bps;
calibRes.RMSE_price   = RMSE_price;
calibRes.RMSE_vol_bps = RMSE_vol;
calibRes.expiry       = expiry_y;
calibRes.tenor        = tenor_y;
calibRes.exitflag     = best_ef;

%% ── 8. Print summary ─────────────────────────────────────────────────────
fprintf('\n--- MHW CALIBRATION RESULTS ---\n');
fprintf('  a     = %8.4f%%   (ref: 12.94%%)\n', a_hat*100);
fprintf('  sigma = %8.4f%%   (ref:  1.26%%)\n', sigma_hat*100);
fprintf('  gamma = %8.5f%%   (ref:  0.07%%)\n', gamma_hat*100);
fprintf('  Price RMSE = %.4e\n', RMSE_price);
fprintf('  Vol   RMSE = %.3f bps\n', RMSE_vol);
fprintf('\n  %-6s %-6s %10s %10s %8s\n','Expiry','Tenor','Mkt(bps)','MHW(bps)','Err(bps)');
fprintf('  %s\n', repmat('-',1,48));
for i = 1:M
    fprintf('  %3dy   %3dy  %8.2f  %8.2f  %+7.3f\n', ...
            expiry_y(i), tenor_y(i), vol_bps(i), vi_bps(i), err_bps(i));
end
fprintf('\n');
end


%% =========================================================================
%%  LOCAL FUNCTION – 
%% =========================================================================
function err2 = swaptionObj(p, expiry_y, tenor_y, OIS, EUR6M, p_mkt)
a=p(1); sig=p(2); gam=p(3);
err2 = 0;
for i = 1:length(expiry_y)
    try
        pv = mhwPrice(a, sig, gam, expiry_y(i), tenor_y(i), OIS, EUR6M);
        if isfinite(pv)
            err2 = err2 + (pv - p_mkt(i))^2;
        else
            err2 = err2 + 1e6;
        end
    catch
        err2 = err2 + 1e6;
    end
end
if ~isfinite(err2),  err2 = length(expiry_y)*1e6;  end
end
