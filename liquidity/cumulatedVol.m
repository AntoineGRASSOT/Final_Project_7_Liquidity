function [Sigma_i, Sigma_N] = cumulatedVol(MHWparams, tau, ti, tN, t0)

a = MHWparams.a;
sigma = MHWparams.sigma;

ti_yf = yearfrac(t0, ti, 3);
tN_yf = yearfrac(t0, tN, 3);
tau_yf = yearfrac(t0, tau, 3);

ksi_i = sigma/a * (1 - exp(-a*(ti_yf - tau_yf)));
ksi_N = sigma/a * (1 - exp(-a*(tN_yf - tau_yf)));

Sigma_i = sqrt(ksi_i.^2 * (1 - exp(-2*a*tau_yf)) / (2*a));
Sigma_N = sqrt(ksi_N^2 * (1 - exp(-2*a*tau_yf)) / (2*a));

end