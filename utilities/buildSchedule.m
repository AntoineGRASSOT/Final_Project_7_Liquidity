function schedule = buildSchedule(t0, te)
% BUILDSCHEDULE  Annual dates from t0 to te included
% The stub (if any) is at the FRONT (short stub in advance, per the paper)

[y0, m0, d0] = datevec(t0);
[ye, ~, ~] = datevec(te);

nYears = ye - y0;    % integer 

schedule = zeros(nYears + 1, 1);
schedule(1) = t0;
for k = 1:nYears
    schedule(k+1) = datenum(y0 + k, m0, d0);
end
% last date is exactly te
schedule(end) = te;
end