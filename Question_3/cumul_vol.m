function cum_vol = cumul_vol(ksi, a, tau)

cum_vol = sqrt((ksi^2) * (1-exp(-2*a*tau))/(2*a));

end