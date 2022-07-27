function Hess = calcHessian(obj)
% calculates the numerical Hessian of function evaluated at a parameter set
% (params) using finite differences.
%
% Adapted from MATLAB code in fminusub, line 431-548.
%
% .........................................................................
% Create: March 9, 2020 (Masaki Tanaka, Bank of Japan)
%
%% --------------------------------------------------------------
% Set the calibration flag to unity (temporarily)
obj.isCalib = 1; 

% Back up the estimated values of the parameters
param_vals_bu = obj.Parameters{1, :};

%% --------------------------------------------------------------
param_vals    = obj.Parameters{1, :};
nParams       = length(param_vals);
f             = obj.runIEKF;
Hess          = zeros(nParams);

% Define stepsize  
Chg = eps^(1/4) * sign(param_vals) .* max(abs(param_vals), 1);

% Make sure step size lies within DiffMinChange and DiffMaxChange
Chg = sign(Chg + eps) .* max(abs(Chg), obj.minDiff);

% Calculate the upper triangle of the finite difference Hessian element 
% by element, using only function values. The forward difference formula 
% we use is
%
% Hessian(i,j) = 1/(h(i)*h(j)) * [f(x+h(i)*ei+h(j)*ej) - f(x+h(i)*ei) 
%                          - f(x+h(j)*ej) + f(x)]                   (2) 
% 
% The 3rd term in (2) is common within each column of Hessian and thus
% can be reused. We first calculate that term for each column and store
% it in the row vector fplus_array.

fplus_array = zeros(1, nParams);

%% --------------------------------------------------------------
for i = 1 : nParams
    xplus                = param_vals;
    xplus(i)             = param_vals(i) + Chg(i);
    obj.Parameters{1, :} = xplus;

    fplus                = obj.runIEKF;       
    fplus_array(i)       = fplus;
end

for i = 1 : nParams
    % For each row, calculate the 2nd term in (4). This term is common to
    % the whole row and thus it can be reused within the current row: we
    % store it in fplus_i.
    xplus                = param_vals;
    xplus(i)             = param_vals(i) + Chg(i);
    obj.Parameters{1, :} = xplus;
    fplus_i              = obj.runIEKF;       

    for j = i : nParams   % start from i: only upper triangle
      % Calculate the 1st term in (2); this term is unique for each element
      % of Hessian and thus it cannot be reused.
        xplus                = param_vals;
        xplus(i)             = param_vals(i) + Chg(i);
        xplus(j)             = xplus(j)      + Chg(j);
        obj.Parameters{1, :} = xplus;
        fplus                = obj.runIEKF;
        
        Hess(i, j) = ( fplus - fplus_i - fplus_array(j) + f )/...
                     ( Chg(i) * Chg(j) ); 
    end 
end

% Fill in the lower triangle of the Hessian
Hess = Hess + triu(Hess, 1)';

%% --------------------------------------------------------------
% Reset the calbiration flag
obj.isCalib = 0; 

% Recover the estimated values of the parameters
obj.Parameters{1, :} = param_vals_bu;
end

