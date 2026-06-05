function [df_OIS, dates_OIS] = bootstrapOIS(dates, data)
% BOOTSTRAPOIS  Bootstrap the OIS discounting curve P^D(t0, ti)
%
% Short end (<=1y): deposit formula
% Long end  (>1y) : OIS swap recursive formula
% ACT/360 year fractions 
%
% Inputs:
%   dates   - what we get with readExcelData_OIS
%             dates.settlement : settlement date
%             dates.ois(i)     : end dates of OIS instruments
%   data    - from readExcelData_OIS
%             data.ois_rates   : [Bid | Ask | Mid], mid = column 3
%
% Outputs:
%   df_OIS    - discount factors P^D(t0, ti), (18x1)
%   dates_OIS - corresponding end dates (same as dates.ois)

t0  = dates.settlement;          
endDates = dates.ois;     % 18 pillar dates
R = data.ois_rates(:, 3);     % Mid quotes (column 3)
N = length(R);

df_OIS    = zeros(N, 1);
dates_OIS = endDates;

for i = 1:N

    te = endDates(i);
    R_i = R(i);
    delta_i = yearfrac(t0, te, 2);   

    if delta_i <= 1.0       
        % P^D(t0, te) = 1 / (1 + delta * R)
        df_OIS(i) = 1 / (1 + delta_i * R_i);

    else
        % Periods: t0 -> t1 -> t2 -> ... -> ti = te included
        schedule = buildSchedule(t0, te);
        n = length(schedule) - 1;    % number of annual periods

        % R * sum_{k=1}^{i-1} delta_k * P^D(t0, t_k)
        annuitySum = 0;
        for k = 1 : n-1
            tk_prev = schedule(k);
            tk      = schedule(k+1);
            delta_k = yearfrac(tk_prev, tk, 2);   % ACT/360

            % Interpolate P^D at tk from already-bootstrapped pillars
            df_k = interpolateOIS(t0, dates_OIS(1:i-1), df_OIS(1:i-1), tk);
            annuitySum = annuitySum + delta_k * df_k;
        end

        delta_n = yearfrac(schedule(end-1), schedule(end), 2);
        df_OIS(i) = (1 - R_i * annuitySum) / (1 + delta_n * R_i);
    end
    
end

fprintf('OIS bootstrap complete: %d pillars\n', N);
fprintf('  df(1w)  = %.8f   z(1w)  = %.4f%%\n', ...
    df_OIS(1), -log(df_OIS(1))/yearfrac(t0,endDates(1),3)*100);
fprintf('  df(12y) = %.8f   z(12y) = %.4f%%\n', ...
    df_OIS(end), -log(df_OIS(end))/yearfrac(t0,endDates(end),3)*100);
end