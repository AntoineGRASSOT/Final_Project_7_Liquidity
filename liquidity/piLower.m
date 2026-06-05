function pi_L = piLower(Sigma_i, Sigma_N)
% this function computes by closed formula the term pi_lower needed in the
% formula for the illiquid bonds
%
% the 2 inputs are the cumulated vol needed for the computation

A = @(eta) exp(-eta/2 .* Sigma_i .* (Sigma_i-Sigma_N));
B = @(eta) 1 + sqrt(pi/2 .* (1-eta)) .* Sigma_N .* exp((1-eta)/8 .* Sigma_N.^2) .*...
    normcdf(sqrt(1-eta)/2 .* Sigma_N);
C = @(eta) 1 + sqrt(pi.*eta/2) .* (2*Sigma_i-Sigma_N) .*...
    exp(eta/8 .* (2*Sigma_i-Sigma_N).^2) .*...
    normcdf(sqrt(eta)/2 .* (2*Sigma_i-Sigma_N));

f = @(eta) 1./(sqrt(1-eta).*sqrt(eta)) .* A(eta) .* B(eta) .* C(eta);

pi_L = exp(-Sigma_N^2/8)/pi .* integral(f, 0, 1, 'RelTol',1e-10, 'AbsTol',1e-12);

end
