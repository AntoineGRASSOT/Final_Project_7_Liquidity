function d = addMonths(t0, k)
% ADDMONTHS  Add k months to date t0 with modified-following convention
%
% Input:
%   t0  - base date (serial)
%   k   - number of months to add (scalar or vector)
%
% Output:
%   d   - resulting date(s) after mod-foll adjustment

k = k(:);
d = zeros(size(k));

[y0, m0, day0] = datevec(t0);

for i = 1:length(k)
    total_months = m0 + k(i);
    y_new = y0 + floor((total_months - 1) / 12);
    m_new = mod(total_months - 1, 12) + 1;

    % End-of-month adjustment: if day0 > days in new month, use last day
    last_day = eomday(y_new, m_new);
    d_new    = min(day0, last_day);

    raw = datenum(y_new, m_new, d_new);
    d(i) = modFoll(raw);
end
end