function pi_U = piUpper(Sigma_i)
% this function computes by closed formula the term pi_upper needed in the
% formula for the illiquid bonds
%
% the inputs are the cumulated vol needed for the computation

pi_U = (4 + Sigma_i.^2)/2 .* normcdf(Sigma_i/2) + ...
    Sigma_i/sqrt(2*pi) .* exp(-Sigma_i.^2/8);

end