function Z = interpolateZSpread(t0, knownDates, knownZ, targetDates)
% Linear on Zeta Spread Z(t) = -ln(Bbar/B)/(t-t0), with B the risk-free OIS  
% discount factor and Bbar the defaultable discount factor

targetDates = targetDates(:);
knownDates = knownDates(:);
knownZ = knownZ(:);

t_known = yearfrac(t0, knownDates, 3);   % ACT/365
t_target = yearfrac(t0, targetDates, 3);

M = length(t_target);
Z = zeros(M, 1);

for k = 1:M
    t = t_target(k);

    if t <= t_known(1)
        % Before first pillar: flat at first Zeta Spread
        Z(k) = knownZ(1);
        continue;
    end

    if t >= t_known(end)
        % After last pillar: flat at last Zeta Spread
        Z(k) = knownZ(end);
        continue;
    end

    Z(k) = interp1(t_known, knownZ, t, "linear");
end

end