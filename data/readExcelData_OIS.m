function [dates, data] = readExcelData_OIS(filename, formatData)
fprintf('--- Loading file : %s ---\n', filename);

% 1.OIS (reference Dates  & Curves OIS)

fprintf('Loading sheet : OIS...\n');

% gloabl Dates (H3 et H5)
tmp_today = readcell(filename, 'Sheet', 'OIS', 'Range', 'H3');
dates.today = smartConvert(tmp_today{1}, formatData);

tmp_settlement = readcell(filename, 'Sheet', 'OIS', 'Range', 'H5');
dates.settlement = smartConvert(tmp_settlement{1}, formatData);

% Dates of the OIS curve (B2:B19)
tmp_ois_dates = readcell(filename, 'Sheet', 'OIS', 'Range', 'B2:B19');
for i = 1:length(tmp_ois_dates)
    dates.ois(i,1) = smartConvert(tmp_ois_dates{i}, formatData);
end

% Numercial data OIS : Bid, Ask, Mid (C2:E19)
% Divide by 100 to change % into decimal (ex: -0.132% -> -0.00132)
data.ois_rates = readmatrix(filename, 'Sheet', 'OIS', 'Range', 'C2:E19') / 100;

data.ois_tenors = readcell(filename, 'Sheet', 'OIS', 'Range', 'A2:A19');


% 2.  Liborvs6m

fprintf('Loading Sheet : Liborvs6m...\n');

% Dates Euribor 6m (B2:B14)
tmp_libor_dates = readcell(filename, 'Sheet', 'Liborvs6m', 'Range', 'B2:B14');
for i = 1:length(tmp_libor_dates)
    dates.libor6m(i,1) = smartConvert(tmp_libor_dates{i}, formatData);
end

% price Euribor 6m (C2:C14)
data.libor6m_prices = readmatrix(filename, 'Sheet', 'Liborvs6m', 'Range', 'C2:C14') / 100;

data.libor6m_tenors = readcell(filename, 'Sheet', 'Liborvs6m', 'Range', 'A2:A14');

% FRA

fprintf('Loading Sheet : FRA...\n');

% start Dates (C2:C10) and end Dates (D2:D10)
tmp_fra_start = readcell(filename, 'Sheet', 'FRA', 'Range', 'C2:C10');
tmp_fra_end   = readcell(filename, 'Sheet', 'FRA', 'Range', 'D2:D10');

for r = 1:size(tmp_fra_start, 1)
    dates.fra(r,1) = smartConvert(tmp_fra_start{r,1}, formatData); % Start
    dates.fra(r,2) = smartConvert(tmp_fra_end{r,1}, formatData);   % End
end

% numerical data FRA : Bid, Ask, Mid (E2:G10)
data.fra_rates = readmatrix(filename, 'Sheet', 'FRA', 'Range', 'E2:G10') / 100;

data.fra_labels = readcell(filename, 'Sheet', 'FRA', 'Range', 'A2:A10');


% Vols

fprintf('Loading Sheet : Vols...\n');

% Vol in bps (D2:D11) 
data.vols = readmatrix(filename, 'Sheet', 'Vols', 'Range', 'D2:D11');


% little check
fprintf('Today       : %s\n', datestr(dates.today));
fprintf('Settlement  : %s\n', datestr(dates.settlement));
fprintf('--- Success in the loading ---\n');
end

function dnum = smartConvert(val, fmt)
    if ismissing(val)
        dnum = NaN;
        return;
    end
    
    if isnumeric(val)
        % Correction du décalage Excel standard (39497 -> 19-Feb-2008)
        dnum = val + datenum('30-Dec-1899');
    elseif isdatetime(val)
        dnum = datenum(val);
    else
        % Si la date est une chaîne de caractères (ex: '21-Sep-15')
        try
            dnum = datenum(val, fmt);
        catch
            try
                dnum = datenum(val);
            catch
                dnum = NaN;
            end
        end
    end
end