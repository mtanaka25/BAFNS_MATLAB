function [R, varargout] = calcSpotRates(obj, factor, ELB, T2M)
% calculates spot rates from the factors and the level of ELB.
% Since spot rates are calculated by averaging forward rates,
% "calcForwardRates.m" is required.
%
% .........................................................................
%  Create March 9, 2020 (Masaki Tanaka)
%

%% 1. Calculate spot rates
Tau2take = round(T2M / obj.dTau); 

% Calculate forward rates
[FR, dFR_dx, info] = obj.calcForwardRates(factor, ELB, T2M);

if info == 1
    R     = NaN;
    dR_dx = NaN;
    varargout{1} = dR_dx;
    varargout{2} = info;
    return
end

% Calculate spot rates by averaging the forward rates
R     = cumsum(FR)     ./ (1 : length(FR))';
dR_dx = cumsum(dFR_dx) ./ repmat((1:length(dFR_dx))', 1, 3);

% Extract the maturities you needed
R     = R(Tau2take);
dR_dx = dR_dx(Tau2take,:);

%% 2. Rearrange optional outputs
varargout{1} = dR_dx;
varargout{2} = info;
end

