function survP = survivalProb(ZS, t0, tau)

tau_yf = yearfrac(t0, tau, 3);
ZS_tau = interpolateZSpread(t0, ZS.T, ZS.Z, tau);
survP = exp(-ZS_tau * tau_yf);

end
