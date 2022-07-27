function setProperties(obj)
% set default properties
%
%.........................................................................
% Create: March 9. 2020 (Masaki Tanaka, Bank of Japan)
%

%% 1. Default Options

obj.nFactors     = 3;
obj.isELBmodel   = 1;
obj.ModelFreq    = 'weekly';

obj.DataFreq = 'daily';

obj.isCalib     = 0;
obj.nIterEKF    = -1e-5;
obj.dTau        = 0.01;
obj.MaxIter     = 100000;
obj.isComputeSE = 1;
obj.minDiff     = 1e-5;

obj.DataFile   = 'JP_data_NZY.csv'; % Name of the data file
obj.T2M        = [0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20];
obj.T2M2Use    = [0.5, 1, 2, 3, 5, 7, 10, 20];
obj.ParamsFile = 'JP_params.csv';

switch obj.ModelFreq
    case 'daily'
        obj.dt = 1/252; % dt equals to one business day
    case 'weekly'
        obj.dt = 1/52;  % dt equals to one week
    case 'monthly'
        obj.dt = 1/12;  % dt equals to one month
    case 'quarterly'
        obj.dt = 1/4;   % dt equals to one quarter
    case 'annually'
        obj.dt = 1;     % dt equals to one year
end

