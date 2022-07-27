function [FR, varargout] = calcForwardRates(obj, factor, ELB, T2M)
% caluculates foward rates from three factors and the level of ELB.
%
% ........................................................................
% Create: March 9, 2020 (Masaki Tanaka, Bank of Japan)
%

%% 0. Loading the required values from the class properties
phi     = obj.Parameters.phi(1)  ;
sig1    = obj.Parameters.sig1(1) ;
sig2    = obj.Parameters.sig2(1) ;
sig3    = obj.Parameters.sig3(1) ;
rho12   = obj.Parameters.rho12(1);
rho13   = obj.Parameters.rho13(1);
rho23   = obj.Parameters.rho23(1);

TauMax  = max(T2M);
Tau_vec = (0 : obj.dTau : TauMax)';

info = 0;
%% 1. Expected short-term rates (under the Q measure)
%  --- In the B-ANSM class, specification of expected short-term rates
%      corresponds to the dynamic arbitrage-free Nelson Siegel model.
 
coef1 = ones(length(Tau_vec), 1); % Coef on the 1st factor
coef2 = exp(- phi * Tau_vec)    ; % Coef on the 2nd factor
coef3 = phi * Tau_vec .* coef2  ; % Coef on the 3rd factor

E_SR_vec = coef1 * factor(1) + coef2 * factor(2) + coef3 * factor(3);

%% 2. Volatility effect (Jensen inequality term)
Omega11 = sig1^2;
Omega22 = sig2^2;
Omega33 = sig3^2;
Omega12 = rho12 * sig1 * sig2;
Omega13 = rho13 * sig1 * sig3;
Omega23 = rho23 * sig2 * sig3;

G1Phi = (1 - coef2) / phi;
F1Phi = G1Phi - Tau_vec .* coef2;

% For the detail of the volatility effect term, see Krippner's textbook
% (Krippner 2014, pp. 77, Equation 3.86)
VE    =  0.5 * Omega11  * Tau_vec .* Tau_vec...
       + 0.5 * Omega22  * G1Phi   .* G1Phi...
       + 0.5 * Omega33 .* F1Phi   .* F1Phi...
       +       Omega12  * Tau_vec .* G1Phi...
       +       Omega13  * Tau_vec .* F1Phi...
       +       Omega23  * G1Phi   .* F1Phi;
   
if any(VE < 0)
    info    = 1;
    FR      = NaN;
    dFR_dx  = NaN;
    varargout{1} = dFR_dx;
    varargout{2} = info;
    return
end

%% 3. Calculating foward rates
% Foward Rates = Expected short-term rates - Jensen inequality
FR = E_SR_vec - VE;

% --------------- [apply the effective lower bound constraint] ------------
if obj.isELBmodel == 1
    g2Phi = coef2 .* coef2;
    G2Phi = 0.5 * (1 - g2Phi) / phi;
    F2Phi = G2Phi - Tau_vec .* g2Phi;
    OmegaSqrd =  Omega11 * Tau_vec...
               + Omega22 * G2Phi...
               + Omega33 * 0.5 * (F2Phi - phi * Tau_vec .* Tau_vec .* g2Phi)...
               + Omega12 * 2   * G1Phi...
               + Omega13 * 2   * F1Phi...
               + Omega23 * F2Phi;
    if any(OmegaSqrd < 0)
        info    = 1;
        FR     = NaN;
        dFR_dx = NaN;
        return
    end
    Omega = sqrt(OmegaSqrd);
    
    d = (FR - ELB) ./ Omega;
    if any(imag(d) ~= 0)
        d(imag(d) ~= 0) = 0;
    end
    normcdf_d = normcdf(d);

    % Calculate gradiant
    dFR_dx = [coef1 .* normcdf_d, ...
              coef2 .* normcdf_d, ...
              coef3 .* normcdf_d];
    % Calculate the ELB-constrained forward rates
    % (See Krippner's textbook (pp. 111-112) for the detail).
    FR    = ELB + (FR - ELB) .* normcdf_d ...
            + exp(-0.5 * d .* d) .* Omega / sqrt(2 * pi);    
        
else
    % If the model is ordinarily affine, the derivatives correspond to 
    % the coefficients in the arbitrage-free Nelson Siegel model.
    dFR_dx = [coef1, coef2, coef3];
end

%%
varargout{1} = dFR_dx;
varargout{2} = info;

end

