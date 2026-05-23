function ksi = compute_ksi(sigma, a, ti, tau)

% CHECK
if ti <= tau
    warning('compute_ksi: ti <= tau detected (ti = %.6f, tau = %.6f)', ti, tau);
end

ksi = (sigma/a)*(1-exp(-a*(ti-tau))) ;

end