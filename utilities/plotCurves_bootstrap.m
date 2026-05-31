function plotCurves_bootstrap(curves)
% PLOTCURVES  Plot OIS and Euribor 6m zero rate curves
%
% Shows:
%   - Zero rates (continuously compounded, ACT/365) for both curves
%   - Discount factors for both curves
%   - Spread: Euribor 6m zero rate minus OIS zero rate

t0 = curves.settle;

% Fine grid for smooth plot (0.1y steps up to 12y)
t_grid  = (0.1 : 0.1 : 12)';
d_grid  = t0 + t_grid * 365;   % approximate serial dates

% Interpolate both curves on fine grid
z_OIS = zeros(size(t_grid));
z_E6m = zeros(size(t_grid));

for k = 1:length(t_grid)
    df_ois = interpolateOIS(t0, curves.OIS.dates, curves.OIS.df, d_grid(k));
    z_OIS(k) = -log(df_ois) / t_grid(k) * 100;

    df_e6m = interpolateOIS(t0, curves.E6m.dates, curves.E6m.df, d_grid(k));
    z_E6m(k) = -log(df_e6m) / t_grid(k) * 100;
end

spread_bps = (z_E6m - z_OIS) * 100;   % in basis points

figure('Name', 'Multicurve Bootstrap - EUR 10-Sep-2015', ...
    'Position', [100 100 1200 400]);

%% Panel 1: Zero rates
subplot(1,3,1);
hold on;
plot(t_grid, z_OIS, 'b-',  'LineWidth', 2);
plot(t_grid, z_E6m, 'r--', 'LineWidth', 2);
% Mark pillar points
plot(curves.OIS.t, curves.OIS.zeroRate*100, 'bs', ...
    'MarkerSize', 6, 'MarkerFaceColor', 'b');
plot(curves.E6m.t, curves.E6m.zeroRate*100, 'r^', ...
    'MarkerSize', 5, 'MarkerFaceColor', 'r');
hold off;
xlabel('Maturity (years)');
ylabel('Zero Rate (%)');
title('Zero Rates (ACT/365, continuous)');
legend('OIS', 'Euribor 6m', 'Location', 'NorthWest');
grid on;
xlim([0 12]);

%% Panel 2: Discount factors
subplot(1,3,2);
hold on;
df_OIS_grid = exp(-z_OIS/100 .* t_grid);
df_E6m_grid = exp(-z_E6m/100 .* t_grid);
plot(t_grid, df_OIS_grid, 'b-',  'LineWidth', 2);
plot(t_grid, df_E6m_grid, 'r--', 'LineWidth', 2);
hold off;
xlabel('Maturity (years)');
ylabel('Discount Factor');
title('Discount / Pseudo-Discount Factors');
legend('P^D(t_0,T) OIS', 'P(t_0,T) Euribor 6m', 'Location', 'NorthEast');
grid on;
xlim([0 12]);

%% Panel 3: Spread E6m - OIS
subplot(1,3,3);
plot(t_grid, spread_bps, 'k-', 'LineWidth', 2);
xlabel('Maturity (years)');
ylabel('Spread (bps)');
title('Euribor 6m - OIS Spread');
grid on;
xlim([0 12]);
yline(0, 'k--');

sgtitle(sprintf('EUR Multicurve Bootstrap - Settlement %s', ...
    datestr(t0, 'dd-mmm-yyyy')));

end