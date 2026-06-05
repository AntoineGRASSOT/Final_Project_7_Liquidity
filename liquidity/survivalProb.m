function survP = survivalProb(ZS, t0, tau)
% SURVIVALPROB Compute survival probability from the Z-spread curve.

tau_yf = yearfrac(t0, tau, 3);
ZS_tau = interpolateZSpread(t0, ZS.T, ZS.Z, tau);
survP = exp(-ZS_tau * tau_yf);

end
