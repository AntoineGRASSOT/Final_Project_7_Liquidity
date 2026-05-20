function d = modFoll(d_input)
% MODFOLL  Apply modified-following business day convention
%
% If d_input is a weekend or holiday, move to next business day.
% If that crosses into the next month, move backward instead.
%

% Input/Output: serial date numbers

d = d_input(:);

for k = 1:length(d)
    dw = weekday(d(k));   % 1=Sun, 7=Sat

    if dw == 7        % Saturday -> move forward 2
        d_cand = d(k) + 2;
    elseif dw == 1    % Sunday   -> move forward 1
        d_cand = d(k) + 1;
    else
        continue;     % already a weekday
    end

    % Check if we crossed into a new month
    [~, m_orig] = datevec(d(k));
    [~, m_new]  = datevec(d_cand);

    if m_new ~= m_orig
        % Modified: go backward to Friday instead
        if dw == 7
            d(k) = d(k) - 1;  % Friday
        else
            d(k) = d(k) - 2;  % Friday
        end
    else
        d(k) = d_cand;
    end
end
end