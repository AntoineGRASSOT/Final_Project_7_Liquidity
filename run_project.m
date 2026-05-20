%% 
clc; clear all;
addpath("data\")
addpath("utilities\")

%%
formatDate = 'dd/mm/yyyy';
% bootstrap datas to compute zero rates
[dates, data] = readExcelData_OIS('curves20150910_project.xlsx', formatDate);
[df_OIS, dates_OIS] = bootstrapOIS(dates, data);

%%
[E6m_dates, E6m_df] = bootstrapEuribor6m(dates, data, dates_OIS, df_OIS);