function runEstimation(obj)
%  Running this script triggers (i) estimating a shadow rate term-structure 
%  model with a time-varying effective lower bound and (ii) implementing
%  simulations using the estimated model.
%
%........................................................................
% Create: March 3, 2020 (Masaki Tanaka, Bank of Japan) 
%


%% **************************************************************
%% 1. Estimate parameters and state variables using the IEKF
t_est_begin = datetime; 

% ------------------------- [1-1 Run Estimation] -------------------------
if obj.isCalib
   %***********************************************************************
   %                      Parameters are calibrated
   %***********************************************************************
    fprintf('*******************************************\n')
    fprintf('  Run the iterated extended Kalman filter  \n')
    fprintf('  using calibrated parameters              \n')
    fprintf('*******************************************\n')
    [minusLLH, factor_mat] = obj.runIEKF;
else
   %***********************************************************************
   %                      Parameters are estimated
   %***********************************************************************
    fprintf('*******************************************\n')
    fprintf('  Run the iterated extended Kalman filter  \n')
    fprintf('  and estimate parameter values            \n')
    fprintf('*******************************************\n')
    [minusLLH, factor_mat] = obj.runMLE; 
end

obj.LLH             = - minusLLH;
obj.FilteredFactors = factor_mat';

% ---------------------- [1-2 Calc standard errors] -----------------------
% If needed, compute standard errors of the estimates, based on the Hessian
if obj.isComputeSE && ~obj.isCalib
    fprintf('%s:: Computing the Hessian... \n', mfilename);
    
    Hessian = obj.calcHessian;
    SE_vec  = sqrt( abs( diag(inv(Hessian)) ) )';
    
    obj.Parameters{2, :} = SE_vec;
end

%% **************************************************************
%% 2. Calculate model-implied spot rate curve
fprintf('%s:: Computing the fitted yield curve... \n', mfilename);
nObs     = size(obj.DateLabel, 1);
FittedYC = NaN(nObs, size(obj.T2M2Use, 2));

for i = 1 : nObs
    R_i = obj.calcSpotRates(obj.FilteredFactors(i, :), obj.Paraeters.ELB, obj.T2M2Use);
    FittedYC(i, :) = R_i;
end

obj.FittedYC = FittedYC ;

%% **************************************************************
obj.Time2Est = datetime - t_est_begin;
end