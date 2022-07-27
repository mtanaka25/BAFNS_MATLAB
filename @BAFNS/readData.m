function readData(obj)
% loads input data 
%
%..........................................................................
% Create: March 3rd, 2020 (Masaki Tanaka, Bank of Japan)
%


%% --------------------------------------------------------------
% load the data file as a "table"
T = readtable(obj.DataFile);

% take date label from the table
if isdatetime(T.DATE)
    obj.DateLabel = T.DATE;
else
    obj.DateLabel = datetime(T.DATE, 'ConvertFrom','datenum');
end

T = removevars(T,{'DATE'});

% take yield curve data from the table
obj.Data = T{:, ismember(obj.T2M, obj.T2M2Use)};

obj.FirstDay = obj.DateLabel(1)  ; % First date of the data
obj.LastDay  = obj.DateLabel(end); % Last date of the data

if isempty(obj.DataFreq)
    % Check the frequency of the original data
    freq_ave = nanmean(obj.DateLabel(2:end) - obj.DateLabel(1:end-1));
    if freq_ave < 7
        obj.DataFreq = 'daily';
    elseif freq_ave < 30
        obj.DataFreq = 'weekly';
    elseif freq_ave < 90
        obj.DataFreq = 'monthly';
    elseif freq_ave < 365
        obj.DataFreq = 'quarterly';
    else
        obj.DataFreq = 'annualy';
    end
end
end

