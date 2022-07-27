function WriteResult2Excel(obj, varargin)
%
%..........................................................................
% Create: March  9, 2020 (Masaki Tanaka, Bank of Japan)
%
%

%% 1. Check the Input ********************************************
fname = sprintf('%s_%s_result.xlsx', datestr(datetime, 'yymmdd'), class(obj));

if nargin == 1
    fprintf(['%s:: You did not declare the name of the output Excel file.',...
        ' The name is automatically set.\n'], mfilename, fname);
else
    if ischar(varargin{1})
        fname = varargin{1};
    else
        fprintf(['%s:: You put an inappropriate argument into %s method. \n',...
                 '                    %s method recieves only one argument and the argument\n'... 
                 '                    should be a character vector of the Excel file to be made.\n'],...
           mfilename, mfilename, mfilename);
    end
end

%% 2. Write out the Result to an Excel File *********************
fprintf('%s:: Saving the estimation and simulation results as an Excel file (%s).\n',...
    mfilename, fname);

% 2-1. Observed and Filtted Yields ---------------------------------------
%
label    = obj.DateLabel;
nSample  = size(label, 1);
zeroline = zeros(nSample, 1);

% Prepare the variable names
varnames = cell(1, size(obj.T2M2Use, 2));
years    = floor(obj.T2M2Use * 12./ 12);  
months   =   rem(obj.T2M2Use * 12 , 12);  
for i = 1 : size(obj.T2M2Use, 2)
    if years(i) == 0
        varnames{1, i} = sprintf('M%d', months(i));
    else
        if months(i) == 0
            varnames{1, i} = sprintf('Y%d', years(i));
        else
            varnames{1, i} = sprintf('Y%dM%d', years(i), months(i));
        end
    end
end 

% Table for Observed Values
Table1_1 = array2table(obj.Data * 100);
Table1_1.Properties.VariableNames = varnames;

% Table for Fitted Values
Table1_2 = array2table(obj.FittedYC * 100);
Table1_2.Properties.VariableNames = varnames;


% Write out
Sheet1   = 'YieldData'; 
Table1_0 = table(label, zeroline,...
    'VariableNames', {'Date', 'Zeroline'});

writetable(Table1_0, fname, 'Sheet', Sheet1,...
    'Range', mtExcelAddress(2, 1),...
    'WriteVariableNames', true);
writetable(Table1_1, fname, 'Sheet', Sheet1,...
    'Range', mtExcelAddress(2, 4),...
    'WriteVariableNames', true);
writetable(Table1_2, fname, 'Sheet', Sheet1,...
    'Range', mtExcelAddress(2, 4 + size(obj.T2M2Use, 2)),...
    'WriteVariableNames', true);

column_label         = cell(1,  size(obj.T2M2Use, 2) + 1);
column_label{1,  1 } = 'Observed_Yields';
column_label{1, end} = 'Fitted_Yields'  ;
xlswrite(fname, column_label, Sheet1, mtExcelAddress(1, 4));

% 2-2. Parameters ---------------------------------------------------------
Sheet2 = 'Parameters';
Table2 = table(obj.Parameters{1, :}',obj.Parameters{2, :}', obj.Parameters{3, :}',...
    'VariableNames', {'PointEstimate', 'SE', 'InitialValue'},...
    'RowNames', obj.Parameters.Properties.VariableNames');

write(Table2, fname, 'Sheet', Sheet2, ...
    'WriteVariableNames', true, 'WriteRowNames', true);
end