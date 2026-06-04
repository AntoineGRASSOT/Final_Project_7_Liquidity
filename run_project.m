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
results = computePremiumDifference(bonds, OIS, ZS_all, MHWparams, issuers, t0, taus, tauLabels);
plotPremiumDifference(bonds, results, issuers, tauLabels);

%% iv) BOND YIELDS (liquid and illiquid 2w / 2m)
yields = computeBondYields(bonds, OIS, ZS_all, issuers, t0, taus, tauLabels, results);
plotBondYields(bonds, yields, issuers, tauLabels);       
plotLiquiditySpread(bonds, yields, issuers, tauLabels);  