%% 
clc; clear all;
addpath("data\")
addpath("utilities\")
addpath("liquidity\")
addpath("calibration\")

%% i) MULTI-CURVE BOOTSRAP (OIS + Euribor 6m)
formatDate = 'dd/mm/yyyy';
[dates, data] = readExcelData_OIS('curves20150910_project.xlsx', formatDate);
[df_OIS, dates_OIS] = bootstrapOIS(dates, data);

[E6m_dates, E6m_df] = bootstrapEuribor6m(dates, data, dates_OIS, df_OIS);

t0 = dates.settlement;  curves.settle = t0;
curves.OIS.dates = dates_OIS;      curves.OIS.df  = df_OIS;
t_OIS = yearfrac(t0, dates_OIS, 3);    % ACT/365
curves.OIS.t = t_OIS;
curves.OIS.zeroRate = -log(df_OIS) ./ t_OIS;
curves.E6m.dates = E6m_dates;    curves.E6m.df = E6m_df;
t_E6m = yearfrac(t0, E6m_dates, 3);
curves.E6m.t = t_E6m;
curves.E6m.zeroRate = -log(E6m_df) ./ t_E6m;
plotCurves_bootstrap(curves);

%% ii) VOLATILITY PARAMETERS' CALIBRATION
[a_hat, sigma_hat, gamma_hat, calibRes] = calibrateMHW(data, dates, dates_OIS, df_OIS, E6m_dates, E6m_df);
MHWparams.a     = a_hat;
MHWparams.sigma = sigma_hat;
MHWparams.gamma = gamma_hat;

%% iii) SHEER LIQUIDITY PREMIUM
% Load corporate bond data (BNPP, Santander)
t0 = dates.settlement ; 
issuers = {'BNPP', 'Santander'};
bonds = loadBondData(t0);
tau_2w = modFoll(addtodate(t0, 14, 'day')); % 2 weeks
tau_2m = modFoll(addtodate(t0, 2, 'month')); % 2 months

taus    = {tau_2w, tau_2m};
tauLabels = {'2 weeks', '2 months'};

%% Z-spread bootstrap per issuer
OIS.T = dates_OIS;
OIS.DF = df_OIS;

ZS_BNPP = bootstrapZSpread(bonds.BNPP, OIS, t0, 'BNPP');
ZS_Santander = bootstrapZSpread(bonds.Santander, OIS, t0, 'Santander');
ZS_all  = {ZS_BNPP, ZS_Santander};  

%% Difference U-L of sheer liquidity premium 
fprintf('\n=== Difference U-L of sheer liquidity premium ===\n');

% Pre-allocate storage for the illiquid bond yields
results = struct();
for i = 1:length(issuers)
    N = length(bonds.(issuers{i}).maturity);
    results.(issuers{i}).DU   = zeros(N, length(taus));
    results.(issuers{i}).DL   = zeros(N, length(taus));
    results.(issuers{i}).diff = zeros(N, length(taus));
end

for i = 1:length(issuers)
    name    = issuers{i};
    bonds_i = bonds.(name);
    ZS_i    = ZS_all{i};
    N       = length(bonds_i.maturity);

    fprintf('\n--- %s ---\n', name);
    fprintf('  %-14s', 'Maturity');
    for j = 1:length(taus)
        fprintf('   %20s', sprintf('diff @ %s', tauLabels{j}));
    end
    fprintf('\n');

    for k = 1:N
        bond_k = struct('couponDates', bonds_i.couponDates{k}, ...
                        'cashflows',   bonds_i.cashflows{k});
        fprintf('  %-14s', datestr(bonds_i.maturity(k), 'dd-mmm-yyyy'));
        for j = 1:length(taus)
            DU = liquidityPremium(bond_k, OIS, ZS_i, MHWparams, t0, taus{j}, 'U');
            DL = liquidityPremium(bond_k, OIS, ZS_i, MHWparams, t0, taus{j}, 'L');

            results.(name).DU(k, j)   = DU;
            results.(name).DL(k, j)   = DL;
            results.(name).diff(k, j) = DU - DL;

            fprintf('   %20.6e', DU - DL);
        end
        fprintf('\n');
    end
end

plotPremiumDifference(bonds, results, issuers, tauLabels);

%% iv) BOND YIELDS (liquid and illiquid 2w / 2m)
% For each bond: liquid yield Y from the dirty model price P, illiquid
% yield y_s from P^s = P - Delta^s (upper bound DU), and liquidity spread
% L_s = (y_s - Y) in basis points
fprintf('\n=== Bond yields (liquid / illiquid) ===\n');

yields = struct();
for i = 1:length(issuers)
    N = length(bonds.(issuers{i}).maturity);
    yields.(issuers{i}).Yliq   = zeros(N, 1);                 % liquid yield
    yields.(issuers{i}).Yill   = zeros(N, length(taus));      % illiquid yields
    yields.(issuers{i}).spread = zeros(N, length(taus));      % L_s in bp
end

for i = 1:length(issuers)
    name    = issuers{i};
    bonds_i = bonds.(name);
    ZS_i    = ZS_all{i};
    N       = length(bonds_i.maturity);

    fprintf('\n--- %s ---\n', name);
    fprintf('  %-14s   %10s', 'Maturity', 'Y (%)');
    for j = 1:length(taus)
        fprintf('   %14s', sprintf('L @ %s (bp)', tauLabels{j}));
    end
    fprintf('\n');

    for k = 1:N
        cd_k = bonds_i.couponDates{k};
        cf_k = bonds_i.cashflows{k};
        bond_k = struct('couponDates', cd_k, 'cashflows', cf_k);

        % Liquid price and liquid yield (full cash-flow set)
        P = bondDirtyModel(bond_k, OIS, ZS_i, t0);
        Y = bondYield(P, cd_k, cf_k, t0);
        yields.(name).Yliq(k) = Y;

        fprintf('  %-14s   %10.4f', datestr(bonds_i.maturity(k), 'dd-mmm-yyyy'), Y * 1e2);

        for j = 1:length(taus)
            Ps = P - results.(name).DU(k, j);      % Illiquid price using the upper-bound premium
            ys = bondYield(Ps, cd_k, cf_k, t0);    % Illiquid yield on the SAME full cash-flow set
            Ls = (ys - Y) * 1e4;                   % Liquidity spread in basis points

            yields.(name).Yill(k, j)   = ys;
            yields.(name).spread(k, j) = Ls;

            fprintf('   %14.4f', Ls);
        end
        fprintf('\n');
    end
end

plotBondYields(bonds, yields, issuers, tauLabels);
