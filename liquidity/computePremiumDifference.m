function results = computePremiumDifference(bonds, OIS, ZS_all, MHWparams, issuers, t0, taus, tauLabels)
% COMPUTEPREMIUMDIFFERENCE  Upper/lower-bound sheer liquidity premium and
% their difference for every bond and every time-to-liquidate (point iii).
%
% Input:
%   bonds      struct with one field per issuer (output of loadBondData)
%   OIS        risk-free curve with fields .T (pillar dates), .DF
%   ZS_all     cell array of Z-spread curves, one per issuer
%   MHWparams  vol params with .a and .sigma
%   issuers    cell array of issuer names (e.g. {'BNPP','Santander'})
%   t0         settlement date (datenum)
%   taus       cell array of time-to-liquidate dates (datenums)
%   tauLabels  cell array of TTL labels (e.g. {'2 weeks','2 months'})
%
% Output:
%   results    struct with one field per issuer, each containing matrices
%              .DU, .DL, .diff (N_bonds x N_ttl) per unit face value

fprintf('\n=== Difference U-L of sheer liquidity premium ===\n');

results = struct();
for i = 1:length(issuers)
    name = issuers{i};
    bonds_i = bonds.(name);
    ZS_i = ZS_all{i};
    N = length(bonds_i.maturity);

    results.(name).DU = zeros(N, length(taus));
    results.(name).DL = zeros(N, length(taus));
    results.(name).diff = zeros(N, length(taus));

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

            results.(name).DU(k, j) = DU;
            results.(name).DL(k, j) = DL;
            results.(name).diff(k, j) = DU - DL;

            fprintf('   %20.6e', DU - DL);
        end
        fprintf('\n');
    end
end

end
