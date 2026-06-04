function yields = computeBondYields(bonds, OIS, ZS_all, issuers, t0, taus, tauLabels, results)
% COMPUTEBONDYIELDS  Liquid and illiquid bond yields and the liquidity yield
% spread for every bond and every time-to-liquidate (point iv).
%
% Input:
%   bonds      struct with one field per issuer (output of loadBondData)
%   OIS        risk-free curve with fields .T (pillar dates), .DF
%   ZS_all     cell array of Z-spread curves, one per issuer
%   issuers    cell array of issuer names (e.g. {'BNPP','Santander'})
%   t0         settlement date (datenum)
%   taus       cell array of time-to-liquidate dates (datenums)
%   tauLabels  cell array of TTL labels (e.g. {'2 weeks','2 months'})
%   results    output of computePremiumDifference (provides .DU)
%
% Output:
%   yields     struct with one field per issuer, each containing
%              .Yliq  (N x 1)        liquid yield (decimal, ACT/365)
%              .Yill  (N x N_ttl)    illiquid yields
%              .spread(N x N_ttl)    liquidity yield spread L_s in bp

fprintf('\n=== Bond yields (liquid / illiquid) ===\n');

yields = struct();
for i = 1:length(issuers)
    name = issuers{i};
    bonds_i = bonds.(name);
    ZS_i = ZS_all{i};
    N = length(bonds_i.maturity);

    yields.(name).Yliq = zeros(N, 1);
    yields.(name).Yill = zeros(N, length(taus));
    yields.(name).spread = zeros(N, length(taus));

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
            Ps = P - results.(name).DU(k, j);   % Illiquid price (upper-bound premium)
            ys = bondYield(Ps, cd_k, cf_k, t0); % Illiquid yield
            Ls = (ys - Y) * 1e4;                % Liquidity yield spread in basis points

            yields.(name).Yill(k, j) = ys;
            yields.(name).spread(k, j) = Ls;

            fprintf('   %14.4f', Ls);
        end
        fprintf('\n');
    end
end

end