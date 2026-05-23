%% 
clc; clear all;
addpath("data\")
addpath("utilities\")
addpath("Question_3\")

%% i) Multi-curve bootstrap (OIS + Euribor 6m)
formatDate = 'dd/mm/yyyy';
[dates, data] = readExcelData_OIS('curves20150910_project.xlsx', formatDate);
[df_OIS, dates_OIS] = bootstrapOIS(dates, data);

[E6m_dates, E6m_df] = bootstrapEuribor6m(dates, data, dates_OIS, df_OIS);

%% ii)
% Calibration volatility parameters
a_hat = ; 
sigma_hat = ;

%% iii) Load corporate bond data (BNPP, Santander)
today = dates.today ; t0 = dates.settlement ; 
bonds = loadBondData(t0);
tau_2w = modFoll(addtodate(today, 14, 'day')); % 2 weeks
tau_2m = modFoll(addtodate(today, 2, 'month')); % 2 months

taus    = {tau_2w, tau_2m};
tauLabels = {'2 weeks', '2 months'};

%% Z-spread bootstrap per issuer
OIS.T = dates_OIS;
OIS.DF = df_OIS;

ZS_BNPP = bootstrapZSpread(bonds.BNPP, OIS, t0, 'BNPP');
ZS_Santander = bootstrapZSpread(bonds.Santander, OIS, t0, 'Santander');




%  Default probabilities from Table 4 
% Values given as 1 - P(t0, tau)
defProb.BNPP.w2  = 1.27e-4;
defProb.BNPP.m2  = 5.42e-4;
defProb.SANT.w2  = 1.80e-4;
defProb.SANT.m2  = 7.73e-4;
survProb.BNPP.w2 = 1 - defProb.BNPP.w2;
survProb.BNPP.m2 = 1 - defProb.BNPP.m2;
survProb.SANT.w2 = 1 - defProb.SANT.w2;
survProb.SANT.m2 = 1 - defProb.SANT.m2;

results = struct();
 
for iIssuer = 1:2
    issuer = issuers{iIssuer};
    % Select Z-spread curve and bond data
    if strcmp(issuer, 'BNPP')
        Z_dates   = Z_dates_BNPP;
        Z_spreads = Z_spreads_BNPP;
        bonds     = bondData.BNPP;
    else
        Z_dates   = Z_dates_SANT;
        Z_spreads = Z_spreads_SANT;
        bonds     = bondData.Santander;
    end
 
    fprintf('\n========================================\n');
    fprintf('  Issuer: %s\n', issuer);
    fprintf('========================================\n');
 
    for iTau = 1:2
        tau      = taus{iTau};
        tauLabel = tauLabels{iTau}; 
        % Survival probability P(t0, tau)
        if strcmp(issuer, 'BNPP')
            if iTau == 1
                P_surv = survProb.BNPP.w2;
            else
                P_surv = survProb.BNPP.m2;
            end
        else
            if iTau == 1
                P_surv = survProb.SANT.w2;
            else
                P_surv = survProb.SANT.m2;
            end
        end
 
        fprintf('\n  tau = %s,  P(t0,tau) = %.8f\n', tauLabel, P_surv);
 
        %  Loop over all bonds of this issuer 
        Delta_lower_total = 0;
        Delta_upper_total = 0;
 
        fprintf('  %-14s  %6s  %8s  %8s  %8s  %10s  %10s\n', ...
            'Maturity','Coupon','B(t0,ti)','pi_L','pi_U', ...
            'contrib_L','contrib_U');
        fprintf('  %s\n', repmat('-', 1, 80));
 
        for b = 1:length(bonds.maturities)
            T_b = bonds.maturities(b);
            K_b = bonds.coupons(b);
 
            % Build payment schedule for this bond
            payDates = buildPaymentDates(t0, T_b); % stub first
            N_b = length(payDates);
 
            % Cash flows
            prevDate = t0;
            c = zeros(N_b, 1);
            for i = 1:N_b
                dcf  = yearfrac(prevDate, payDates(i), 1);  % Act/Act -> paper
                c(i) = K_b * dcf * 100;
                prevDate = payDates(i);
            end
            c(end) = c(end) + 100;
 
            % Keep only ti > tau
            mask = payDates > tau;
            payDates_filt = payDates(mask);
            c_filt = c(mask);
            N_filt = sum(mask);
 
            if N_filt == 0
                fprintf('  Bond %s: no cash flows after tau, skipped\n', ...
                    datestr(T_b,'dd-mmm-yyyy'));
                continue;
            end
 
            tau_yf = yearfrac(t0, tau, 3);
            % sum of contributions
            for i = 1:N_filt
                ti   = payDates_filt(i); ti_yf = yearfrac(t0,ti,3);
                tN   = payDates_filt(end);  tN_yf = yearfrac(t0,tN,3);
                ci   = c_filt(i);
 
                % B(t0, ti): risk-free ois df
                B_ti = interpolateOIS(t0, dates_OIS, df_OIS, ti);
 
                % pi_i^L and pi_i^U
                pi_low = compute_pi_lower(ti_yf,tN_yf,tau_yf,sigma_hat, a_hat);
                pi_up  = compute_pi_upper(ti_yf,tau_yf,sigma_hat, a_hat);
                
                contrib_L = ci * B_ti * (pi_low - P_surv);
                contrib_U = ci * B_ti * (pi_up - P_surv);
 
                Delta_lower_total = Delta_lower_total + contrib_L;
                Delta_upper_total = Delta_upper_total + contrib_U;
 
                % Print only for first bond as example
                if b == 1
                    fprintf('  %-14s  %6.4f  %8.6f  %8.6f  %8.6f  %10.6f  %10.6f\n', ...
                        datestr(ti,'dd-mmm-yyyy'), ci, B_ti, piL, piU, ...
                        contrib_L, contrib_U);
                end
            end
        end
 
        % Results
        spread_diff = Delta_upper_total - Delta_lower_total;
 
        fprintf('\n  --- Results: %s, tau = %s ---\n', issuer, tauLabel);
        fprintf('  Lower bound Delta_tau^L = %.6f  (%.4f bps)\n', ...
            Delta_lower_total, Delta_lower_total * 1e4);
        fprintf('  Upper bound Delta_tau^U = %.6f  (%.4f bps)\n', ...
            Delta_upper_total, Delta_upper_total * 1e4);
        fprintf('  Difference (U - L)      = %.6f  (%.4f bps)\n', ...
            spread_diff, spread_diff * 1e4);
 
        fn_tau = strrep(tauLabel, ' ', '_');
        results.(issuer).(fn_tau).lower = Delta_lower_total;
        results.(issuer).(fn_tau).upper = Delta_upper_total;
        results.(issuer).(fn_tau).diff  = spread_diff;
    end
end
 
%  Summary table
fprintf('\n\n  ============ SUMMARY TABLE ============\n');
fprintf('  %-12s  %-10s  %12s  %12s  %12s\n', ...
    'Issuer', 'tau', 'Lower (bps)', 'Upper (bps)', 'Diff (bps)');
fprintf('  %s\n', repmat('-', 1, 65));
 
for iIssuer = 1:2
    issuer = issuers{iIssuer};
    for iTau = 1:2
        fn_tau = strrep(tauLabels{iTau}, ' ', '_');
        r = results.(issuer).(fn_tau);
        fprintf('  %-12s  %-10s  %12.4f  %12.4f  %12.4f\n', ...
            issuer, tauLabels{iTau}, ...
            r.lower*1e4, r.upper*1e4, r.diff*1e4);
    end
end

plotLiquidityBounds(results, issuers, tauLabels);
