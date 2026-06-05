function payDates = buildPaymentDates(t0, maturity)
% stub different from build schedule like 15 set 2015 -> 27 nov 2015 -> 27 nov2016 -> 27 nov 2017 
% backward from the maturity date
%
% INPUTS
%   t0        : valuation date
%   maturity  : final maturity date of the instrument
%
% OUTPUT
%   payDates  : column vector containing all annual payment dates

[yM, mM, dM] = datevec(maturity);
[y0, ~,  ~ ] = datevec(t0);
nYears = yM - y0;
payDates = zeros(nYears, 1);
for k = 1:nYears
    raw = datenum(yM - (nYears - k), mM, dM);
    payDates(k) = modFoll(raw);
end
payDates = payDates(payDates > t0);
if ~isempty(payDates)
    payDates(end) = maturity;
end
end
