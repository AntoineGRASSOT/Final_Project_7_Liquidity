function pi_U = piUpper(Sigma_i)

pi_U = (4 + Sigma_i.^2)/2 .* normcdf(Sigma_i/2) + ...
    Sigma_i/sqrt(2*pi) .* exp(-Sigma_i.^2/8);

end