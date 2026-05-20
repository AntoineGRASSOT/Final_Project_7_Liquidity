function df = interpolateOIS(t0, knownDates, knownDFs, targetDates)
% Linear on zero rates z(t) = -ln(P)/(t-t0), ACT/365

targetDates = targetDates(:);
knownDates  = knownDates(:);
knownDFs    = knownDFs(:);
 
t_known  = yearfrac(t0, knownDates, 3);   % ACT/365
t_target = yearfrac(t0, targetDates, 3);
 
% Zero rates at pillars: z = -ln(P) / t
z_known = -log(knownDFs) ./ t_known;
 
M  = length(t_target);
df = zeros(M, 1);
 
for k = 1:M
    t = t_target(k);
 
    if t <= 0
        df(k) = 1.0;
        continue;
    end
 
    if t <= t_known(1)
        % Before first pillar: flat at first zero rate
        df(k) = exp(-z_known(1) * t);
        continue;
    end
 
    if t >= t_known(end)
        % After last pillar: flat at last zero rate
        df(k) = exp(-z_known(end) * t);
        continue;
    end

    z = interp1(t_known, z_known, t, 'linear');
    df(k) = exp(-z * t);
end

end