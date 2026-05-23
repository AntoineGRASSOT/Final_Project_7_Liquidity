function [a_hat, sigma_hat, gamma_hat, calibRes] = calibrateMHW(data, dates, dates_OIS, df_OIS, E6m_dates, E6m_df)
%% CALIBRATEMHW  MHW calibration – 9 diagonal ATM CS swaptions, Baviera (2019).

%% 0. Build curve structs
t0       = dates.settlement;
OIS.T    = yearfrac(t0, dates_OIS(:), 3);
OIS.DF   = df_OIS(:);
EUR6M.T  = yearfrac(t0, E6m_dates(:), 3);
EUR6M.PD = E6m_df(:);

fprintf('\n>>> MHW calibration: 9 diagonal ATM swaptions\n');
fprintf('   OIS  curve: %d pillars, T in [%.4f, %.4f]\n', ...
        numel(OIS.T), min(OIS.T), max(OIS.T));
fprintf('   EUR6M curve: %d pillars, T in [%.4f, %.4f]\n', ...
        numel(EUR6M.T), min(EUR6M.T), max(EUR6M.T));

%% 1. Swaption grid
M        = 9;
%expiry_y = (1:M)';
%tenor_y  = (M:-1:1)';
%%%%%%%%%
vols_raw = data.vols(:);
first_ok = find(isfinite(vols_raw), 1, 'first');   % index of first finite vol

if isempty(first_ok)
    error('calibrateMHW: data.vols contains no finite values. Check the Excel file.');
end
if numel(vols_raw) < first_ok + M - 1
    error('calibrateMHW: not enough finite vols in data.vols (need %d, found %d).', ...
          M, numel(vols_raw) - first_ok + 1);
end

vol_mkt_bps = vols_raw(first_ok : first_ok + M - 1);   % (9x1) bps
vol_mkt     = vol_mkt_bps / 1e4;                        % decimal

expiry_y = (1:M)';    % 1y, 2y, ..., 9y
tenor_y  = (M:-1:1)'; % 9y, 8y, ..., 1y

fprintf('   Vol range read: data.vols(%d:%d)\n', first_ok, first_ok+M-1);
fprintf('   Vols (bps): ');
fprintf('%.2f  ', vol_mkt_bps);
fprintf('\n');
%%%%%%%


vol_mkt  = data.vols(1:M) / 1e4;

%% 2. Market prices
price_mkt = zeros(M,1);
S_atm_vec = zeros(M,1);
B_alpha_v = zeros(M,1);
C_S0_vec  = zeros(M,1);
for i = 1:M
    B_a          = getDF(OIS, expiry_y(i));
    [S0,~,~]     = swapRateATM(expiry_y(i), tenor_y(i), OIS, EUR6M);
    C0           = cashAnnuity(S0, tenor_y(i), 1);
    price_mkt(i) = bachelierCS(S0, S0, expiry_y(i), B_a, C0, vol_mkt(i));
    B_alpha_v(i) = B_a;  S_atm_vec(i) = S0;  C_S0_vec(i) = C0;
end
fprintf('   Market prices: [%.4e , %.4e]\n', min(price_mkt), max(price_mkt));

%% 3. Step-by-step debug at p0 = [0.10, 0.01, 0.001] ----------------------
fprintf('\n   --- Step-by-step debug at p0=[0.10, 0.01, 0.001] ---\n');
a0=0.10; s0=0.01; g0=0.001;
for i = 1:M
    fprintf('   Swaption %dy%dy: ', expiry_y(i), tenor_y(i));
    try
        [pr, xs, Sa] = mhwPrice(a0, s0, g0, expiry_y(i), tenor_y(i), OIS, EUR6M);
        fprintf('price=%.4e  xstar=%.4f  S_atm=%.4f%%\n', pr, xs, Sa*100);
    catch ME
        fprintf('ERROR -> %s\n', ME.message);
    end
end
fprintf('\n');

%% 4. Pre-check objective at each p0
p0s = [0.10, 0.010, 0.001;
       0.05, 0.015, 0.0005;
       0.20, 0.008, 0.002;
       0.15, 0.012, 0.0];

fprintf('   Pre-check objective:\n');
for k = 1:size(p0s,1)
    val_k = swaptionObjective(p0s(k,:), expiry_y, tenor_y, OIS, EUR6M, price_mkt);
    fprintf('   p0(%d) [%.3f %.4f %.5f] -> obj=%.4e  fin=%d\n', ...
            k, p0s(k,1), p0s(k,2), p0s(k,3), val_k, isfinite(val_k));
end

%% 5. fmincon
lb  = [1e-4, 1e-4, 0.0];
ub  = [1.0,  0.10, 1.0];
opt = optimoptions('fmincon','Algorithm','interior-point','Display','off', ...
    'MaxIterations',1000,'MaxFunctionEvaluations',5000, ...
    'FunctionTolerance',1e-12,'StepTolerance',1e-12,'OptimalityTolerance',1e-12);
objFun   = @(p) swaptionObjective(p, expiry_y, tenor_y, OIS, EUR6M, price_mkt);
best_val = Inf;  best_p = p0s(1,:);  best_ef = -1;
fprintf('\n');
for k = 1:size(p0s,1)
    try
        [pk,fk,efk] = fmincon(objFun,p0s(k,:),[],[],[],[],lb,ub,[],opt);
        fprintf('   Start %d: a=%6.4f sig=%6.4f gam=%7.5f Err=%.3e (flag=%d)\n', ...
                k,pk(1),pk(2),pk(3),sqrt(fk),efk);
        if fk<best_val, best_val=fk; best_p=pk; best_ef=efk; end
    catch ME
        fprintf('   Start %d: fmincon error -> %s\n',k,ME.message);
    end
end
a_hat=best_p(1); sigma_hat=best_p(2); gamma_hat=best_p(3);

%% 6. Model prices and implied vols
price_model=zeros(M,1); vol_impl_bps=zeros(M,1);
for i=1:M
    price_model(i)  = mhwPrice(a_hat,sigma_hat,gamma_hat,expiry_y(i),tenor_y(i),OIS,EUR6M);
    vol_impl_bps(i) = bachelierCS_implVol(S_atm_vec(i),S_atm_vec(i),expiry_y(i), ...
                        B_alpha_v(i),C_S0_vec(i),price_model(i))*1e4;
end
vol_mkt_bps=vol_mkt*1e4; error_bps=vol_impl_bps-vol_mkt_bps;
RMSE_price=sqrt(best_val/M); RMSE_vol=sqrt(mean(error_bps.^2));

%% 7. Output
calibRes.price_mkt=price_mkt; calibRes.price_model=price_model;
calibRes.vol_mkt_bps=vol_mkt_bps; calibRes.vol_impl_bps=vol_impl_bps;
calibRes.error_bps=error_bps; calibRes.RMSE_price=RMSE_price;
calibRes.RMSE_vol_bps=RMSE_vol; calibRes.expiry=expiry_y;
calibRes.tenor=tenor_y; calibRes.exitflag=best_ef;

fprintf('\n--- MHW CALIBRATION RESULTS ---\n');
fprintf('  a     = %8.4f%%   (ref: 12.94%%)\n', a_hat*100);
fprintf('  sigma = %8.4f%%   (ref:  1.26%%)\n', sigma_hat*100);
fprintf('  gamma = %8.5f%%   (ref:  0.07%%)\n', gamma_hat*100);
fprintf('  Price RMSE = %.4e\n', RMSE_price);
fprintf('  Vol   RMSE = %.3f bps\n', RMSE_vol);
fprintf('\n  %-6s %-6s %10s %10s %8s\n','Expiry','Tenor','Mkt(bps)','MHW(bps)','Err(bps)');
fprintf('  %s\n',repmat('-',1,48));
for i=1:M
    fprintf('  %3dy   %3dy  %8.2f  %8.2f  %+7.3f\n', ...
            expiry_y(i),tenor_y(i),vol_mkt_bps(i),vol_impl_bps(i),error_bps(i));
end
fprintf('\n');
end

%% =========================================================================
%  LOCAL FUNCTION
%% =========================================================================
function err2 = swaptionObjective(p, expiry_y, tenor_y, OIS, EUR6M, price_mkt)
a=p(1); sig=p(2); gam=p(3);
err2=0;
for i=1:length(expiry_y)
    try
        pv = mhwPrice(a, sig, gam, expiry_y(i), tenor_y(i), OIS, EUR6M);
        if isfinite(pv)
            err2 = err2 + (pv - price_mkt(i))^2;
        else
            err2 = err2 + 1e6;
        end
    catch
        err2 = err2 + 1e6;
    end
end
if ~isfinite(err2), err2=9e6; end
end
