function plotLiquidityBounds(results, issuers, tauLabels)

figure('Name', 'Liquidity Spread Bounds', 'Position', [100 100 900 400]);

colors = {'b', 'r'};
x = 1:2;   % two tau values

for iIssuer = 1:2
    issuer = issuers{iIssuer};
    subplot(1, 2, iIssuer);
    hold on;

    lower_vals = zeros(2,1);
    upper_vals = zeros(2,1);
    for iTau = 1:2
        fn_tau = strrep(tauLabels{iTau}, ' ', '_');
        lower_vals(iTau) = results.(issuer).(fn_tau).lower * 1e4;
        upper_vals(iTau) = results.(issuer).(fn_tau).upper * 1e4;
    end

    bar_data = [lower_vals, upper_vals];
    b = bar(x, bar_data, 'grouped');
    b(1).FaceColor = [0.3 0.5 0.9];
    b(2).FaceColor = [0.9 0.3 0.3];

    set(gca, 'XTickLabel', tauLabels);
    ylabel('\Delta_\tau (bps)');
    title(sprintf('%s - Liquidity Bounds', issuer));
    legend('Lower bound', 'Upper bound', 'Location', 'NorthWest');
    grid on;
    hold off;
end

sgtitle('Sheer Liquidity Premium \Delta_\tau - Upper and Lower Bounds');
end
