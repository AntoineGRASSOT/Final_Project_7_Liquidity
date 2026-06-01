function y = bondYield(price, couponDates, cashflows, t0)
% BONDYIELD  Continuously-compounded yield (IRR) of a coupon bond.
%
% Input:
%   price        target price (per unit face value)
%   couponDates  column of future payment dates 
%   cashflows    column of cash-flows 
%   t0           settlement date
%
% Output:
%   y            continuously-compounded yield (decimal, ACT/365)

tau = yearfrac(t0, couponDates, 3);          % ACT/365

pvFun = @(r) sum(cashflows .* exp(-r .* tau)) - price;

% Bracketed root search (yields expected well inside this interval)
y = fzero(pvFun, [-0.05, 0.30]);

end
