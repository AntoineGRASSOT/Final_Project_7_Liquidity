function plotLiquiditySpread(bonds, yields, issuers, tauLabels)
% PLOTLIQUIDITYSPREAD  Plot the liquidity yield spread L_s (in basis points)
% versus maturity for each issuer, for each time-to-liquidate. 
%
% Input:
%   bonds      struct with one field per issuer (output of loadBondData)
%   yields     struct with one field per issuer, each containing
%              .spread (N_bonds x N_ttl) liquidity yield spread in bp
%   issuers    cell array of issuer names (e.g. {'BNPP','Santander'})
%   tauLabels  cell array of TTL labels (e.g. {'2 weeks','2 months'})

for i = 1:length(issuers)
    name = issuers{i};
    mat = bonds.(name).maturity;

    figure;
    plot(mat, yields.(name).spread(:, 1), 'bs-', ...
         'LineWidth', 1.3, 'MarkerFaceColor', 'b'); hold on;
    plot(mat, yields.(name).spread(:, 2), 'r^--', ...
         'LineWidth', 1.3, 'MarkerFaceColor', 'r');

    yrs = year(min(mat)):(year(max(mat)) + 1);   % include first and last year
    xt = datenum(yrs', 1, 1);                    % one tick per year
    set(gca, 'XTick', xt);
    xlim([xt(1), xt(end)]);
    datetick('x', 'yyyy', 'keepticks', 'keeplimits');
    grid on;
    xlabel('Maturity');
    ylabel('Liquidity yield spread L^s (bp)');
    title(sprintf('Liquidity yield spread - %s', name));
    legend(sprintf('TTL = %s', tauLabels{1}), ...
           sprintf('TTL = %s', tauLabels{2}), ...
           'Location', 'NorthEast');
end

end
