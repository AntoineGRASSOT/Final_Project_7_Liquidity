function plotBondYields(bonds, yields, issuers, tauLabels)
% PLOTBONDYIELDS  Plot liquid and illiquid bond yields versus maturity for each issuer
%
% Input:
%   bonds      struct with one field per issuer (output of loadBondData)
%   yields     struct with one field per issuer, each containing
%              .Yliq  (N_bonds x 1)        liquid yield
%              .Yill  (N_bonds x N_ttl)    illiquid yields (one col per TTL)
%   issuers    cell array of issuer names (e.g. {'BNPP','Santander'})
%   tauLabels  cell array of TTL labels (e.g. {'2 weeks','2 months'})

for i = 1:length(issuers)
    name = issuers{i};
    mat  = bonds.(name).maturity;

    figure;
    plot(mat, yields.(name).Yliq * 1e2, 'ko-', ...
         'LineWidth', 1.3, 'MarkerFaceColor', 'k'); hold on;
    plot(mat, yields.(name).Yill(:, 1) * 1e2, 'bs-', ...
         'LineWidth', 1.3, 'MarkerFaceColor', 'b');
    plot(mat, yields.(name).Yill(:, 2) * 1e2, 'r^--', ...
         'LineWidth', 1.3, 'MarkerFaceColor', 'r');
    yrs = year(min(mat)):(year(max(mat)) + 1);   % include first and last year
    xt  = datenum(yrs', 1, 1);                    % one tick per year
    set(gca, 'XTick', xt);
    xlim([xt(1), xt(end)]);
    datetick('x', 'yyyy', 'keepticks', 'keeplimits');
    grid on;
    xlabel('Maturity');
    ylabel('Yield (%)');
    title(sprintf('Liquid and illiquid bond yields - %s', name));
    legend('Liquid', ...
           sprintf('Illiquid, TTL = %s', tauLabels{1}), ...
           sprintf('Illiquid, TTL = %s', tauLabels{2}), ...
           'Location', 'NorthWest');
end

end