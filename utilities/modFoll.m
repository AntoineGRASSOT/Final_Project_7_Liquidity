function d = modFoll(d_input)
% MODFOLL  Apply modified-following business day convention
%          with TARGET calendar (weekends + ECB/TARGET holidays).
%
% TARGET fixed holidays (every year):
%   01-Jan  New Year's Day
%   01-May  Labour Day
%   25-Dec  Christmas Day
%   26-Dec  Boxing Day
%
% TARGET moveable holidays (computed via Anonymous Gregorian algorithm):
%   Good Friday   (Friday before Easter)
%   Easter Monday (Monday after Easter)
%
% If the adjusted date crosses into the next calendar month,
% the modified-following convention moves BACKWARD to the last
% good business day in the original month.
%
% Input/Output: MATLAB serial datenum (scalar or column vector)

d = d_input(:);

for k = 1:length(d)
    % Advance until we land on a business day
    while ~isBusinessDay(d(k))
        d(k) = d(k) + 1;
    end

    % Check if we crossed into a new month vs. original date
    [~, m_orig] = datevec(d_input(k));
    [~, m_new ] = datevec(d(k));

    if m_new ~= m_orig
        % Modified-following: step backward until we find a business day
        d(k) = d_input(k) - 1;
        while ~isBusinessDay(d(k))
            d(k) = d(k) - 1;
        end
    end
end
end


%% =========================================================================
function ok = isBusinessDay(dn)
% Returns true if dn (serial datenum) is neither a weekend nor a TARGET holiday.

dw = weekday(dn);
if dw == 1 || dw == 7   % Sunday or Saturday
    ok = false;
    return;
end

[y, m, day] = datevec(dn);

% Fixed TARGET holidays
fixed = [datenum(y,  1,  1);   % New Year's Day
         datenum(y,  5,  1);   % Labour Day
         datenum(y, 12, 25);   % Christmas
         datenum(y, 12, 26)];  % Boxing Day

if any(dn == fixed)
    ok = false;
    return;
end

% Moveable TARGET holidays: Good Friday and Easter Monday
[gf, em] = easterDates(y);
if dn == gf || dn == em
    ok = false;
    return;
end

ok = true;
end


%% =========================================================================
function [goodFriday, easterMonday] = easterDates(y)
% Anonymous Gregorian algorithm for Easter Sunday.

a = mod(y, 19);
b = floor(y / 100);
c = mod(y, 100);
d = floor(b / 4);
e = mod(b, 4);
f = floor((b + 8) / 25);
g = floor((b - f + 1) / 3);
h = mod(19*a + b - d - g + 15, 30);
i = floor(c / 4);
k = mod(c, 4);
l = mod(32 + 2*e + 2*i - h - k, 7);
m = floor((a + 11*h + 22*l) / 451);
month  = floor((h + l - 7*m + 114) / 31);
day_e  = mod(h + l - 7*m + 114, 31) + 1;

easter       = datenum(y, month, day_e);
goodFriday   = easter - 2;
easterMonday = easter + 1;
end
