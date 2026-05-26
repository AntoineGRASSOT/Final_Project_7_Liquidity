function plotPremiumDifference(bonds, results, issuers, tauLabels)
% PLOTPREMIUMDIFFERENCE  Plot the difference between upper- and lower-bound
% sheer liquidity premium for each issuer, reproducing Fig. 2 (BNPP) and
% Fig. 3 (Santander) of Baviera-Nassigh-Nastasi (2021).
%
% Input:
%   bonds      struct with one field per issuer (output of loadBondData)
%   results    struct with one field per issuer, each containing a matrix
%              .diff (N_bonds x N_ttl) of U-L differences in face value
%   issuers    cell array of issuer names (e.g. {'BNPP','Santander'})
%   tauLabels  cell array of TTL labels (e.g. {'2w','2m'})

for i = 1:length(issuers)
    name = issuers{i};
    figure;
    plot(bonds.(name).maturity, results.(name).diff(:, 1), 'bs-', ...
         'LineWidth', 1.3, 'MarkerFaceColor', 'b'); hold on;
    plot(bonds.(name).maturity, results.(name).diff(:, 2), 'r^--', ...
         'LineWidth', 1.3, 'MarkerFaceColor', 'r');
    datetick('x', 'yyyy'); grid on;
    xlabel('Maturity');
    ylabel('\Delta^U - \Delta^L  (face value)');
    title(sprintf('Upper-lower bound difference of sheer liquidity premium - %s', name));
    legend(sprintf('TTL = %s', tauLabels{1}), ...
           sprintf('TTL = %s', tauLabels{2}), ...
           'Location', 'NorthWest');
end

end
