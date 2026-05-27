function C = cashAnnuity(S, n, m)
%% CASHANNUITY  Cash annuity C_{alpha,omega}(S) for CS swaptions.
%
%    C(S) = sum_{i=1}^{n*m} (1/m)/(1+S/m)^i
%         = (1/S) * (1 - 1/(1+S/m)^(n*m)),  S > -m.
%
%  Limit for S -> 0 is n. The fixed leg is annual, so the default is m=1.

if nargin < 3
    m = 1;
end

N      = n * m;
orig_sz = size(S);          
S      = S(:);           

singular = abs(S) < 1e-10;
invalid  = S <= -m & ~singular;
regular  = ~singular & ~invalid;

C = NaN(size(S));
C(singular) = N / m;
C(regular)  = (1 ./ S(regular)) .* (1 - 1 ./ (1 + S(regular)/m).^N);
% C(invalid) stays NaN

C = reshape(C, orig_sz);   % restore the row

end
