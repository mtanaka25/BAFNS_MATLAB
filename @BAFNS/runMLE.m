function [minusLLH, x_T] = runMLE(obj)
% estimates Black-type arbitrage-free Nelson-Siegel (B-AFNS) model with 
% three factors by the method of maximum likelihood estimation.
%
%  * For the detail of the framework, see Krippner [2013].
%
% This code is called as a method of BAFNS class.
%
% [required methods]
%  * runIEKF
%    --- implements the iterated extended Kalman filter.
%..........................................................................
% Create: January 23, 2020 (Masaki Tanaka, Bank of Japan)
% Update: March    2, 2020 (Masaki Tanaka, Bank of Japan)
% --- Enable to impose a restriction of [-1, 1] on three "rho"s
%    more efficiently
% Update: March    6, 2020 (Masaki Tanaka, Bank of Japan)
% --- Make the whole system object-oriented
%

%% 1. Prepare for the optimization 
% convert rho12, rho13, and rho23 in order to keep their values between
% the interval of (-1, 1).
obj.Parameters.rho12 = obj.Parameters.rho12 ./(1 - obj.Parameters.rho12);
obj.Parameters.rho13 = obj.Parameters.rho13 ./(1 - obj.Parameters.rho13);
obj.Parameters.rho23 = obj.Parameters.rho23 ./(1 - obj.Parameters.rho23);

% Move the initial values to the third row of the parameter table 
obj.Parameters{3, :} = obj.Parameters{1, :}; 

% Set options for optimization loutines. 
Options    = optimset('Display'    , 'iter'      , 'TolFun' , 1e-2       ,...
                      'TolX'       , Inf         , 'MaxIter', obj.MaxIter,...
                      'MaxFunEvals', obj.MaxIter);
 
%% 2. Run a MATLAB optimization routine
% Minimize "minus log-likelihood" using MATLAB's fminsearch     
[OptParams, minusLLH, Exitflag, Output] = ...
    fminsearch(@ObjFunc, obj.Parameters{3, :}, Options);

%% 3. Save the results
% Store the estimated parameters in the first row in the parameter table
obj.Parameters{1, :} = OptParams;

% Convert rho12, rho13, and rho23 
obj.Parameters.rho12  = obj.Parameters.rho12 ./(1 + abs(obj.Parameters.rho12));
obj.Parameters.rho13  = obj.Parameters.rho13 ./(1 + abs(obj.Parameters.rho13));
obj.Parameters.rho23  = obj.Parameters.rho23 ./(1 + abs(obj.Parameters.rho23));

% Store the stats of the estimation
Output.exitflag   = Exitflag;
obj.EstStats      = Output;

% Calc filtered factors using the estimated parameters
obj.isCalib = 1; % Turn on the calibration mode temporarily
[~, x_T]           = obj.runIEKF;
obj.isCalib = 0; % Reset the estimation/calibration flag

%% --------------------------------------------------------------
%% Log-likelihood function (nested in the "runMLE" routine)
    function LLH = ObjFunc(param_vals)
        obj.Parameters{1, :} = param_vals;       
        [LLH, ~]           = obj.runIEKF;
    end
end