function [E6m_dates, E6m_df] = bootstrapEuribor6m(dates, data, OIS_dates, OIS_df)
% BOOTSTRAPEURIBOR6M  Bootstrap the Euribor 6m pseudo-discount curve

t0 = dates.settlement;  

% Monthly knots t̂_1 ... t̂_12 (mod-foll)
t_hat = addMonths(t0, (1:12)');   % (12x1)

t_hat18 = addMonths(t0, 18);   % t̂_18 (18 months)

% Annual knots 2y...12y: use the IRS end dates from data
% These are dates.libor6m(2:end) -> 1y is index 1, 2y is index 2, etc.
% We skip 6m (index 1) and 1y (index 2) since those are already handled
IRS_endDates = dates.libor6m(:);     % 13 dates: 6m,1y,2y,...,12y
IRS_rates    = data.libor6m_prices;  % 13 rates (already /100)

% Full set of curve knots (will grow as we bootstrap)
% Initialize with t0 (pseudo-df = 1 by definition)
knot_dates = t0;
knot_df    = 1.0;

% STEP 1: 6m deposit -> P(t̂_6)
% P(t0, t̂_6) = 1 / (1 + delta(t0, t̂_6) * L^6m)
% L^6m = IRS_rates(1) (6m Euribor deposit rate)

L_6m    = IRS_rates(1);    % 6m Euribor rate
t6      = IRS_endDates(1);    % t̂_6 = 14-Mar-2016
delta_6 = yearfrac(t0, t6, 2);   % ACT/360
P_6     = 1.0 / (1.0 + delta_6 * L_6m);

knot_dates = [knot_dates; t6];
knot_df    = [knot_df;    P_6];

fprintf('  Step 1: P(t0, t̂_6)  = %.8f   [6m depo = %.4f%%]\n', P_6, L_6m*100);

% STEP 2: 6×12 FRA -> P(t̂_12)   [Forward move]
% FRA 6x12: start = t̂_6, end = t̂_12
% F_fra = rate of this FRA from data.fra_rates
% P(t0; t̂_6, t̂_12) = 1 / (1 + delta(t̂_6, t̂_12) * F_fra)
% P(t0, t̂_12) = P(t0, t̂_6) * P(t0; t̂_6, t̂_12)

% 6x12 FRA is in our data  
% FRA start dates are in dates.fra(:,1), end dates in dates.fra(:,2)
% The 6x12 FRA starts at ~t̂_6 and ends at ~t̂_12
t12     = IRS_endDates(2);   % 1y end date = t̂_12 = 14-Sep-2016

fra_start = dates.fra(:,1);
fra_end   = dates.fra(:,2);
fra_mid   = data.fra_rates(:,3);   % Mid = column 3

% Match FRA whose start date is closest to t̂_6
[~, idx_6x12] = min(abs(fra_start - t6));
F_6x12    = fra_mid(idx_6x12);
delta_612 = yearfrac(t6, t12, 2);                % ACT/360
P_6to12   = 1.0 / (1.0 + delta_612 * F_6x12);   % pseudo-df [t̂_6, t̂_12]
P_12      = P_6 * P_6to12;

knot_dates = [knot_dates; t12];
knot_df    = [knot_df;    P_12];

fprintf('  Step 2: P(t0, t̂_12) = %.8f   [6x12 FRA = %.4f%%]\n', P_12, F_6x12*100);

% STEP 3: FRAs 1×7 ... 5×11 -> P(t̂_1) ... P(t̂_5)   [Backward move]

% For i = 1, 2, 3, 4, 5:
%   - FRA i×(i+6): F(t0; t̂_i, t̂_{i+6}) from data
%   - Interpolate P(t0, t̂_{i+6}) from {t̂_6, t̂_12} already known
%   - P(t̂_i) = P(t̂_{i+6}) * (1 + delta(t̂_i, t̂_{i+6}) * F_fra)

% Note: we do NOT yet have t̂_1...t̂_5 on the knot grid,
% so interpolation uses the two known knots t̂_6 and t̂_12.

fprintf('  Step 3: Backward move from FRAs 1x7 ... 5x11\n');

for i = 1:5
    t_i    = t_hat(i);       % t̂_i
    t_ip6  = t_hat(i+6);     % t̂_{i+6}

    % Interpolate P(t0, t̂_{i+6}) using current knots
    P_ip6 = interpolateOIS(t0, knot_dates, knot_df, t_ip6);

    % Find matching FRA: start ~= t̂_i, end ~= t̂_{i+6}
    [~, idx_fra] = min(abs(fra_start - t_i));
    F_i = fra_mid(idx_fra);

    delta_i = yearfrac(t_i, t_ip6, 2);     % delta(t̂_i, t̂_{i+6}): ACT/360

    % Backward formula: P(t̂_i) = P(t̂_{i+6}) / P(t0; t̂_i, t̂_{i+6})
    %                            = P(t̂_{i+6}) * (1 + delta_i * F_i)
    P_i = P_ip6 * (1.0 + delta_i * F_i);

    knot_dates = [knot_dates; t_i];
    knot_df    = [knot_df;    P_i];

    fprintf('    i=%d: t̂_%d: P = %.8f  [FRA %dx%d = %.4f%%]\n', ...
        i, i, P_i, i, i+6, F_i*100);
end

% Sort knots by date
[knot_dates, si] = sort(knot_dates);
knot_df = knot_df(si);

% STEP 5: 18m and 2y swaps -> F_3, F_4  (no STIR futures for 6m)

% f = 2 (semi-annual floating leg)
% Forward rates: F_1 = F(t0; t0, t̂_6)    already known
%                F_2 = F(t0; t̂_6, t̂_12)  already known
%                F_3 = F(t0; t̂_12, t̂_18) UNKNOWN -> from 18m swap
%                F_4 = F(t0; t̂_18, t_2y)  UNKNOWN -> from 2y swap

% Floating leg dates (semi-annual, mod-foll):
%   t̃_1 = t̂_6,  t̃_2 = t̂_12,  t̃_3 = t̂_18,  t̃_4 = t_2y

% I(1) = S_18m * [delta(t0,t̂_6)*P^D(t0,t̂_6) + delta(t̂_6,t̂_18)*P^D(t0,t̂_18)]
%   (fixed leg annuity of 18m swap: semi-annual with two periods)

% I(2) from 2y swap, f=2:
%   w_1*F_1 + w_2*F_2 + w_3*F_3 + w_4*F_4 = I(2)
%   => F_4 = (I(2) - I(1)) / w_4

fprintf('  Step 5: 18m and 2y swaps -> F_3, F_4\n');

% ---- Floating leg payment dates (semi-annual) ----
% t̃_1 = t̂_6, t̃_2 = t̂_12, t̃_3 = t̂_18, t̃_4 = t_2y
t_2y   = IRS_endDates(3);   % 14-Sep-2017
t_flt  = [t_hat(6); t12; t_hat18; t_2y];   % 4 payment dates

% OIS discount factors at floating leg payment dates
PD_flt = interpolateOIS(t0, OIS_dates, OIS_df, t_flt);

% delta~ for floating leg (ACT/360)
flt_starts = [t0; t_flt(1:end-1)];
delta_flt  = yearfrac(flt_starts, t_flt, 2);   % ACT/360

% w_k = delta~_k * P^D(t0, t̃_k)
w = delta_flt .* PD_flt;   % (4x1)

% ---- Known forward rates F_1, F_2 ----
% F_1 = F(t0; t0, t̂_6): from pseudo-dfs at t0 and t̂_6
%   P(t0; t0, t̂_6) = P(t0,t̂_6)/P(t0,t0) = P_6
%   F_1 = (1/P_6 - 1) / delta(t0, t̂_6)
delta_01 = yearfrac(t0, t_hat(6), 2);
P_t0     = 1.0;   % pseudo-df at t0 = 1
F_1 = (P_t0 / P_6 - 1.0) / delta_01;

% F_2 = F(t0; t̂_6, t̂_12)
delta_12_fwd = yearfrac(t_hat(6), t12, 2);
F_2 = (P_6 / P_12 - 1.0) / delta_12_fwd;

fprintf('    F_1 = %.6f%%  F_2 = %.6f%%\n', F_1*100, F_2*100);

% ---- 18m swap: I(1) ----
% Fixed leg of 18m swap: semi-annual periods (f=2, so 2 periods for 18m? 
% for 18m swap with f=2, fixed leg has 3 semi-annual periods:
%   [t0,t̂_6], [t̂_6,t̂_12], [t̂_12,t̂_18]
% But standard EUR IRS fixed leg is ANNUAL. 
% "fixed leg is paid annually" -> for 18m, one annual payment at 1y + stub
% Actually for <2y swaps, the fixed leg has the same frequency as floating
% In the 6m case: both fixed and floating are semi-annual for short maturities
%
% I(1) = S_18m * annuity_fixed(18m)
% The fixed leg annuity uses the OIS discount factors:
%   annuity = sum delta_k^fix * P^D(t0, T_k^fix)
%
% For 18m: fixed dates = t̂_6, t̂_12, t̂_18 (semi-annual, same as float for <2y)
% delta_k^fix: ACT/360 for EUR IRS

% 18m swap rate: interpolate between 1y (IRS(2)) and 2y (IRS(3))
% since we don't have 18m explicitly -> linear interpolation on swap rates
t_18m_yrs = yearfrac(t0, t_hat18, 3);
t_1y_yrs  = yearfrac(t0, IRS_endDates(2), 3);
t_2y_yrs  = yearfrac(t0, IRS_endDates(3), 3);
w_interp  = (t_18m_yrs - t_1y_yrs) / (t_2y_yrs - t_1y_yrs);
S_18m     = IRS_rates(2) * (1 - w_interp) + IRS_rates(3) * w_interp;

fprintf('    S_18m (interpolated) = %.6f%%\n', S_18m*100);

% Fixed leg dates for 18m swap (semi-annual)
fix_dates_18m = [t_hat(6); t12; t_hat18];
PD_fix_18m    = interpolateOIS(t0, OIS_dates, OIS_df, fix_dates_18m);
fix_starts_18m = [t0; fix_dates_18m(1:end-1)];
delta_fix_18m  = yearfrac(fix_starts_18m, fix_dates_18m, 2);   % ACT/360
annuity_18m    = sum(delta_fix_18m .* PD_fix_18m);

I1 = S_18m * annuity_18m;
fprintf('    I(1) = %.8f\n', I1);

% Note: w_3 corresponds to t̃_3 = t̂_18
%  we don't have P^D(t̂_18) yet -> interpolate from OIS
PD_t18  = interpolateOIS(t0, OIS_dates, OIS_df, t_hat18);
delta_flt_3 = yearfrac(t12, t_hat18, 2);
w3 = delta_flt_3 * PD_t18;

F_3 = (I1 - w(1)*F_1 - w(2)*F_2) / (2*w3);
fprintf('    F_3 = %.6f%%\n', F_3*100);

% New pseudo-discount at t̂_18:
% P(t0, t̂_18) = P(t0, t̂_12) / (1 + delta(t̂_12, t̂_18) * F_3)
delta_12_18 = yearfrac(t12, t_hat18, 2);
P_18 = P_12 / (1.0 + delta_12_18 * F_3);

knot_dates = [knot_dates; t_hat18];
knot_df    = [knot_df;    P_18];
[knot_dates, si] = sort(knot_dates);
knot_df = knot_df(si);

fprintf('    P(t0, t̂_18) = %.8f\n', P_18);

% ---- 2y swap: I(2) -> F_4 ----
% I(2) = S_2y * sum_{k=1}^{2} delta_k^fix * P^D(t0, T_k^fix)
% Fixed leg of 2y IRS (annual): T_1=t_1y, T_2=t_2y
S_2y = IRS_rates(3);   % 2y IRS rate

fix_dates_2y = [IRS_endDates(2); IRS_endDates(3)];   % 1y, 2y
PD_fix_2y    = interpolateOIS(t0, OIS_dates, OIS_df, fix_dates_2y);
fix_starts_2y = [t0; fix_dates_2y(1)];
delta_fix_2y  = yearfrac(fix_starts_2y, fix_dates_2y, 2);
annuity_2y    = sum(delta_fix_2y .* PD_fix_2y);

I2 = S_2y * annuity_2y;
fprintf('    I(2) = %.8f\n', I2);

% F_4 = (I(2) - I(1)) / w_4
% w_4 = delta(t̂_18, t_2y) * P^D(t0, t_2y)
PD_t2y      = interpolateOIS(t0, OIS_dates, OIS_df, t_2y);
delta_flt_4 = yearfrac(t_hat18, t_2y, 2);
w4 = delta_flt_4 * PD_t2y;

F_4 = (I2 - I1) / w4;
fprintf('    F_4 = %.6f%%\n', F_4*100);

% New pseudo-discount at t_2y:
delta_18_2y = yearfrac(t_hat18, t_2y, 2);
P_2y = P_18 / (1.0 + delta_18_2y * F_4);

knot_dates = [knot_dates; t_2y];
knot_df    = [knot_df;    P_2y];
[knot_dates, si] = sort(knot_dates);
knot_df = knot_df(si);

fprintf('    P(t0, t_2y)  = %.8f\n', P_2y);

% STEP 6: IRS vs 6m, maturities 3y ... 12y

% For each annual maturity t_i (i = 3,...,12):

%   I(i) = S(t0, t_i) * sum_{k=1}^{i} delta_k * P^D(t0, t_k)
%          where delta_k = yearfrac(t_{k-1}, t_k, ACT/360)  [fixed leg, annual]
%          and P^D = OIS discount factors

%   sum_{k=1}^{2i} w_k * F_k = I(i)
%   Last unknown: F_{2i} = (I(i) - sum_{k=1}^{2i-1} w_k*F_k) / w_{2i}

%   Then:
%   P(t0, t_i) = P(t0, t_{i-1}) / [(1 + delta~_{2i-1}*F_{2i-1}) *
%                                    (1 + delta~_{2i}  *F_{2i}  )]
%   where the product covers the two 6m periods of the ith year.

fprintf('  Step 6: IRS 3y ... 12y\n');

% IRS rates for 3y to 12y: 
IRS_long_dates = IRS_endDates(4:end);    % 3y...12y (10 maturities)
IRS_long_rates = IRS_rates(4:end);

% We need to track ALL semi-annual forward rates F_k and their w_k weights
% Collect them in order: F_1...F_4 already known
% F_k corresponds to period [t_{k-1}^flt, t_k^flt] (semi-annual)

% Build the full floating leg schedule up to 12y (semi-annual, mod-foll)
% t̃_0 = t0, t̃_1=t̂_6, t̃_2=t̂_12, t̃_3=t̂_18, t̃_4=t_2y, ..., t̃_24=t_12y
t_flt_full = addMonths(t0, (6:6:144)');   % 6m,12m,...,144m (6 to 144)
% Force the annual knots to match exactly the IRS end dates
for k = 1:length(IRS_endDates)
    % replace t_flt_full entries close to each annual knot
    [minval, idx] = min(abs(t_flt_full - IRS_endDates(k)));
    if minval < 10   % within 10 days
        t_flt_full(idx) = IRS_endDates(k);
    end
end

% All floating leg times (including t0 as start)
t_flt_all = [t0; t_flt_full];

% Forward rates array: F(k) = F(t0; t̃_{k-1}, t̃_k)
nF_total = length(t_flt_full);   % = 24 for 12y
F_all    = zeros(nF_total, 1);
F_all(1) = F_1;
F_all(2) = F_2;
F_all(3) = F_3;
F_all(4) = F_4;

for yr = 3:12    % 3y to 12y
    i = yr;
    t_mat = IRS_long_dates(yr - 2);   % maturity date

    % Fixed leg annual dates and OIS dfs
    fix_schedule = buildSchedule(t0, t_mat);
    PD_fix = interpolateOIS(t0, OIS_dates, OIS_df, fix_schedule(2:end));
    fix_starts_i = fix_schedule(1:end-1);
    delta_fix_i  = yearfrac(fix_starts_i, fix_schedule(2:end), 2);
    annuity_i    = sum(delta_fix_i .* PD_fix);

    I_i = IRS_long_rates(yr-2) * annuity_i;

    % Floating leg weights w_k = delta~_k * P^D(t0, t̃_k) for k=1..2i
    n_flt_periods = 2 * i;   % number of semi-annual periods up to t_i

    % Compute sum of known terms: sum_{k=1}^{2i-1} w_k * F_k
    sumWF = 0.0;
    for k = 1:n_flt_periods - 1
        t_prev_k = t_flt_all(k);
        t_next_k = t_flt_all(k+1);
        delta_k  = yearfrac(t_prev_k, t_next_k, 2);
        PD_k     = interpolateOIS(t0, OIS_dates, OIS_df, t_next_k);
        w_k      = delta_k * PD_k;

        % F_k: if not yet computed, interpolate from current pseudo-dfs
        if k <= 4 || F_all(k) ~= 0
            Fk = F_all(k);
        else
            % Interpolate forward rate from current knots
            P_prev = interpolateOIS(t0, knot_dates, knot_df, t_prev_k);
            P_next = interpolateOIS(t0, knot_dates, knot_df, t_next_k);
            Fk = (P_prev/P_next - 1.0) / delta_k;
            F_all(k) = Fk;
        end
        sumWF = sumWF + w_k * Fk;
    end

    % Last unknown F_{2i}
    k_last   = n_flt_periods;
    t_prev_L = t_flt_all(k_last);
    t_next_L = t_flt_all(k_last + 1);
    delta_L  = yearfrac(t_prev_L, t_next_L, 2);
    PD_L     = interpolateOIS(t0, OIS_dates, OIS_df, t_next_L);
    w_last   = delta_L * PD_L;

    F_last = (I_i - sumWF) / w_last;
    F_all(k_last) = F_last;

    % New pseudo-discount at t_mat:
    % P(t0, t_i) = P(t0, t_{i-1}) * (1/(1+d_{2i-1}*F_{2i-1})) * (1/(1+d_{2i}*F_{2i}))
    P_prev_mat = interpolateOIS(t0, knot_dates, knot_df, t_flt_all(k_last));
    k_pen      = k_last - 1;
    t_prev_pen = t_flt_all(k_pen);
    t_next_pen = t_flt_all(k_pen + 1);
    delta_pen  = yearfrac(t_prev_pen, t_next_pen, 2);
    F_pen      = F_all(k_pen);

    % P(t_{i-1}) already on knot grid
    P_im1 = interpolateOIS(t0, knot_dates, knot_df, t_flt_all(k_last));

    % Two semi-annual steps
    P_new = P_im1 / (1.0 + delta_pen  * F_pen) / (1.0 + delta_L * F_last);

    knot_dates = [knot_dates; t_mat];
    knot_df    = [knot_df;    P_new];
    [knot_dates, si] = sort(knot_dates);
    knot_df = knot_df(si);

    t_mat_yrs = yearfrac(t0, t_mat, 3);
    z_new     = -log(P_new) / t_mat_yrs * 100;
    fprintf('    %2dy: P = %.8f  z = %7.4f%%  [S=%s%%]\n', ...
        yr, P_new, z_new, num2str(IRS_long_rates(yr-2)*100, '%.4f'));
end

mask = knot_dates > t0;  % t0 not useful
E6m_dates = knot_dates(mask);
E6m_df    = knot_df(mask);

fprintf('\n  Euribor 6m Pseudo-Discount Curve:\n');
fprintf('  %-16s  %10s  %10s\n', 'Date', 'P(t0,T)', 'ZeroRate');
fprintf('  %s\n', repmat('-', 1, 42));
for i = 1:length(E6m_dates)
    t_i = yearfrac(t0, E6m_dates(i), 3);
    z_i = -log(E6m_df(i)) / t_i * 100;
    fprintf('  %-16s  %10.7f  %9.4f%%\n', ...
        datestr(E6m_dates(i), 'dd-mmm-yyyy'), E6m_df(i), z_i);
end
end