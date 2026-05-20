function C = cashAnnuity(S, n, m)
%% CASHANNUITY  Cash annuity C_{alpha,omega}(S) for CS swaptions.
%
%  Formula eq. (2.17) in Baviera (2019):
%    C(S) = sum_{i=1}^{n*m} (1/m)/(1+S/m)^i
%         = (1/S) * (1 - 1/(1+S/m)^(n*m)),  S ~= 0.
%
%  For S -> 0 the limit is n. In the EUR market the fixed leg is annual,
%  so the default is m=1.

if nargin < 3
    m = 1;
end

N = n * m;
C = zeros(size(S));

for k = 1:numel(S)
    s = S(k);
    if abs(s) < 1e-10
        C(k) = N / m;
    elseif s <= -m
        C(k) = NaN;
    else
        C(k) = (1 / s) * (1 - 1 / (1 + s / m)^N);
    end
end
end
