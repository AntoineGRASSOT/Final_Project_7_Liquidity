function pi_up = compute_pi_upper(ti,tau, sigma, a)

ksi_i = compute_ksi(sigma,a, ti, tau);
cum_vol = cumul_vol(ksi_i, a ,tau);
pi_up = 0.5 * (4+cum_vol^2) * normcdf(cum_vol/2) + exp((-cum_vol^2)/8) * cum_vol/(sqrt(2*pi)) ;

end