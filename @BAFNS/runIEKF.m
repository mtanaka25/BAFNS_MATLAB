function [minusLLH, x_mat] = runIEKF(obj)
% implements the iterated extended Kalman filter algorithm, which is a
% generalized version of the iterated Kalman filter.  When you set nIterELF
% to unity, this method implements the ordinary extended Kalman filter.
%
%  * For the theoritical detail of these filtering methods, see textbooks
%    about state-space models such as Simon [2006].
%
% This code is called as a method of BAFNS (Black-type argitrage-free 
% Nelson-Siegel) class.
%
% [required methods]
%  * calcSpotRates
%    --- calcurates spot rates in the B-AFNS class.
%  * calcFowardRates
%    --- calcurates instantaneous foward rates in the B-AFNS class.
%    --- This method is called in "calcSpotRates" method.
%
%.........................................................................
% Create: March 6, 2020 (Masaki Tanaka, Bank of Japan)
%

%% Load the required values from the class properties.

data        = obj.Data;
nFactors    = obj.nFactors;
param_names = obj.Parameters.Properties.VariableNames;
dt          = obj.dt;

T2M         = obj.T2M2Use;
isCalib     = obj.isCalib;
nIterEKF    = obj.nIterEKF;
dTau        = obj.dTau;

param_vals  = obj.Parameters{1, :};


%% 1. Set Parameters
[T, K] = size(data); % T: Length of the time-series, K: # of the maturities

% decompose the parameter vector
for i = 1 : length(param_vals)
   eval(sprintf('%s = param_vals(%d);', param_names{i}, i));
end

% transition matrix of the factors under the P measure 
kappaP = [kappaP11, kappaP12, kappaP13;...
          kappaP21, kappaP22, kappaP23;...
          kappaP31, kappaP32, kappaP33];

% steady-state level of the factors 
thetaP = [thetaP1;...
          thetaP2;...
          thetaP3];
      
% size of the factor-specific shocks
sig1 = abs(sig1);
sig2 = abs(sig2);
sig3 = abs(sig3);

if ~isCalib == 1
    % convert rho12, rho13, and rho23 in order to keep their values
    % between the interval of (-1, 1).
    rho12  = rho12/(1 + abs(rho12));
    rho13  = rho13/(1 + abs(rho13));
    rho23  = rho23/(1 + abs(rho23));
end

% 
Sp12      = sqrt(1 - rho12^2);
p123_Sp12 = (rho23 - rho12 * rho13)/Sp12;
SIGMA     = [sig1       , 0             , 0                                  ;...
             rho12*sig2 , sig2*Sp12     , 0                                  ;...
             rho13*sig3 , sig3*p123_Sp12, sig3*sqrt(1-rho13^2-p123_Sp12^2)];
OMEGA     = SIGMA * SIGMA';

% size of measurement errors
sig_nu    = param_vals(end - K + 1 : end);

%% 2. Initializeing Iterated Extended Kalman Filter
% Allocate memory for historical values and variances of the factors
x_mat = NaN(nFactors, T);          
P_mat = NaN(nFactors, nFactors, T); 

% Extended Kalman filter if IEKF_Count=0, iterated EKF with fixed
% IEKF_Count iterations if IEKF_Count>0, and iterated EKF with fixed
% tolerance abs(IEKF_Count) if IEKF_Count<0.

if nIterEKF < 0
    x_Tolerance = abs(nIterEKF);
    isTolerance = 1;
    nIterEKF    = 20;
else
    isTolerance = 0;    
end

% ----------------- [Initializeing state equation] -----------------------
% Calculate the state equation quantities based on parameter values.
[V, D] = eig(kappaP);

% Check if the model is stationary.
d1     = D(1,1);
d2     = D(2,2);
d3     = D(3,3);
if any([real(d1), real(d2), real(d3)] < 0)
    minusLLH   = 10^10;
    fprintf('%s:: The model is not stationary...\n', mfilename);
    return
end

F = expm(-kappaP * dt);

if any(abs(eig(F)) > 1)
    minusLLH   = 10^10;
    return
end

Q = [ G(2*d1, dt), G(d1+d2,dt), G(d1+d3,dt) ;...
      0          , G(2*d2,dt) , G(d2+d3,dt) ;...
      0          , 0          , G(2*d3 ,dt) ];
Q = Q  + transpose(triu(Q, 1));% NOTE: Use TRANSPOSE, because "'" gives conjugate transpose.
U = V  \ OMEGA / transpose(V);
Q = U .* Q;
Q = V  * Q * transpose(V);
Q = real(Q);

% -------------- [Initializeing observation equation] --------------------
% The spec of observation equation depends on the states, thus has to be
% re-calculated at each step of the Kalman filter.

% Matrix of the sizes of measurement errors
R = diag(sig_nu.^2);

% -------------------------- [Initial values] ----------------------------
% x_f and P_f are forecasted values, that is, x(t|t-1) and P(t|t-1)
% x_u and P_u are  updated   values, that is, x(t|t)   and P(t|t)

% The starting values of the factors are their unconditional mean.
x_u = thetaP;

% The initial variances of the factors are their unconditional variances. 
P_u = [0.5/d1, 1/(d1+d2), 1/(d1+d3);...
       0     , 0.5/d2   , 1/(d2+d3);...
       0     , 0        , 0.5/d3   ];
P_u = P_u + transpose(triu(P_u, 1));
P_u = U .* P_u;
P_u = V  * P_u * transpose(V);
P_u = real(P_u);
% Note: Prevent any imaginary parts that arise from round-off errors.

% The log-likelihood starts from zero.
logL = 0;

%% 3. Run Iterated Extended Kalman Filter
for t = 1 : T
    % ---------------------- [Forecasting step] ---------------------------
    x_f =(eye(nFactors) - F) * thetaP + F * x_u;
    P_f = F * P_u * F' + Q;
    
    % ------------------------ [Updating step] ----------------------------
    y_t     = data(t, :)';
    x_u_pre = x_f;
    x_u_now = x_f;
    for i = 1 : 1 + nIterEKF
        % EKF and IEKF iterations
        % y_fitted = h(x_f, 0) , i.e. fitted values of y_t given forecasted x.
        % Jacob    = dR/dX(x_f), i.e. the Jacobian given forecasted x.
        [y_fitted, Jacob_i, iserror] = obj.calcSpotRates(x_u_now, ELB, T2M);
        if iserror
            x_mat    = NaN(nFactors, T);
            minusLLH = 10^10;
            return
        end
        HPHR_i  = (Jacob_i * P_f * Jacob_i' + R);
        K_i     = P_f * Jacob_i' / HPHR_i;
        w_i     = y_t - y_fitted - Jacob_i * (x_f - x_u_now);
        x_u_new = x_f + K_i * w_i;
        if isTolerance
            if all(abs(x_u_new - x_u_now) < x_Tolerance)
                % Difference from last update within tolerance, so exit.
                break
            end
            if all(abs(x_u_new - x_u_pre) < x_Tolerance)
                % Allows for numerical cycling between i+1, i, i-1 updates.
                % Difference from i-1 update within tolerance, so exit.
                x_u_new = 0.5 * (x_u_new + x_u_now);
                break
            end
        end
        % Record these values to allow testing for convergence.
        x_u_pre = x_u_now;
        x_u_now = x_u_new;
    end
    
    % Calculate final posterior values and record values.
    x_u             = x_u_new;
    P_u             = (eye(nFactors) - K_i * Jacob_i) * P_f;
    x_mat(:, t)     = x_u;
    P_mat(:, :, t)  = P_u;
    logL            = logL + log(det(HPHR_i)) + w_i' / HPHR_i * w_i;
end

%% 4. Calculating the log likelihood

% The value to be maximized
LLH = - 0.5 * K * T * log(2 * pi) - 0.5 * logL;

% The value to be "miniimized"
minusLLH = - LLH;

end

%% subfunction
function out = G(phi, tau)
    if phi <= 0
        out = tau;
    else
        out = 1 / phi * (1 - exp(- phi * tau) );
    end
end
