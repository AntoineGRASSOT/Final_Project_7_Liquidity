function bonds = loadBondData(t0)
% LOADBONDDATA  Reference data for BNPP and Santander liquid bonds
%              on 10-Sep-2015 (value date).
% Input:
%   t0   - settlement date 
%
% Output:
%   bonds.<issuer>.maturity   - column of MATLAB datenums
%   bonds.<issuer>.coupon     - column of annual coupon rates (decimal)
%   bonds.<issuer>.cleanPrice - column of clean prices (per unit face)
%   bonds.<issuer>.invoice     - column of invoice prices (clean + accrual)
%

% --- BNPP (Table 2) -----------------------------------------------------
bnpp_mat = datenum({'27-Nov-2017'; ...
                    '12-Mar-2018'; ...
                    '21-Nov-2018'; ...
                    '28-Jan-2019'; ...
                    '23-Aug-2019'; ...
                    '13-Jan-2021'; ...
                    '24-Oct-2022'; ...
                    '20-May-2024'}, 'dd-mmm-yyyy');

bnpp_cpn = [2.875; 1.500; 1.375; 2.000; 2.500; 2.250; 2.875; 2.375] / 100;

bnpp_px  = [105.575; 102.768; 102.555; 104.536; ...
            106.927; 106.083; 110.281; 106.007] / 100;

bonds.BNPP.maturity   = bnpp_mat;
bonds.BNPP.coupon     = bnpp_cpn;
bonds.BNPP.cleanPrice = bnpp_px;

% --- Santander (Table 3) ------------------------------------------------
san_mat = datenum({'27-Mar-2017'; ...
                   '04-Oct-2017'; ...
                   '15-Jan-2018'; ...
                   '20-Apr-2018'; ...
                   '14-Jan-2019'; ...
                   '13-Jan-2020'; ...
                   '24-Jan-2020'; ...
                   '14-Jan-2022'; ...
                   '10-Mar-2025'}, 'dd-mmm-yyyy');

san_cpn = [4.000; 4.125; 1.750; 0.625; 2.000; 0.875; 4.000; 1.125; 1.125] / 100;

san_px  = [105.372; 107.358; 102.766;  99.885; 103.984; ...
            99.500; 112.836;  98.166;  93.261] / 100;

bonds.Santander.maturity   = san_mat;
bonds.Santander.coupon     = san_cpn;
bonds.Santander.cleanPrice = san_px;

% --- Coupon schedules and accruals --------------------------------------
issuers = {'BNPP', 'Santander'};

for i = 1:length(issuers)
    name = issuers{i};
    N    = length(bonds.(name).maturity);

    couponDates = cell(N, 1);
    cashflows   = cell(N, 1);
    accrual     = zeros(N, 1);

    for k = 1:N
        T = bonds.(name).maturity(k);
        c = bonds.(name).coupon(k);

        nmax = ceil((T - t0)/365) + 1;                 % Upper bound 
        candidates = addMonths(T, -12*(0:nmax)');      % Bacward candidtaes dates for the coupons payments
        n_futures = sum(candidates > t0);  
        paymentsDates = modFoll(flipud(candidates(1:n_futures)));
        prevDates = modFoll(candidates(n_futures + 1));

        % Compute the cash flows
        cf = c * ones(length(paymentsDates), 1);
        cf(end) = c + 1;

        % Compute the accruals
        accr = c * ((t0 - prevDates)/(paymentsDates(1) - prevDates));

        couponDates{k} = paymentsDates;
        cashflows{k} = cf;
        accrual(k) = accr;
    end

    bonds.(name).couponDates = couponDates;
    bonds.(name).cashflows = cashflows;
    bonds.(name).accrual = accrual;
    bonds.(name).invoice = bonds.(name).cleanPrice + accrual;
end

% --- Print tables (paper-style, with invoice price column) --------------
fprintf('\nBond data loaded (value date = %s):\n', datestr(t0, 'dd-mmm-yyyy'));

issuerLabels = {'BNPP (Table 2)', 'Santander (Table 3)'};
header_fmt = '  %-14s    %12s    %14s    %16s\n';
row_fmt    = '  %-14s    %12.3f    %14.3f    %16.3f\n';
dash       = repmat('-', 1, 72);

for i = 1:length(issuers)
    name = issuers{i};
    N    = length(bonds.(name).maturity);

    fprintf('\n%s\n', issuerLabels{i});
    fprintf(header_fmt, 'maturity', 'coupon (%)', 'clean price', 'invoice price');
    fprintf('  %s\n', dash);

    for k = 1:N
        fprintf(row_fmt, ...
            datestr(bonds.(name).maturity(k), 'dd-mmm-yyyy'), ...
            bonds.(name).coupon(k)     * 100, ...
            bonds.(name).cleanPrice(k) * 100, ...
            bonds.(name).invoice(k)    * 100);
    end
end

end  