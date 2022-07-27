classdef BAFNS < handle
    % BANSM3 is the MATLAB class for the Black-type arbitrage-free Nelson
    % Siegel model (B-AFNS).
    %
    % ....................................................................
    % Create March 9, 2020 (Masaki Tanaka, Bank of Japan)
    %
    %
    
    %% A. Properties
    properties
        % Mode Spec ******************************************************
        nFactors     % # of factors
        
        FilteredFactors % History of Factors
        
        Parameters      % table of parameters
 
        isELBmodel   % If the model belongs to the shadow-rate term 
                     % structure models, isELBmodel should be set to 1
                       
        ModelFreq    % Frequency of the discretized model
        
        dt           % marginal change of time in the discretized model

        DateLabel    % Date label
                
        FittedYC     % Fitted yield curve
      
        % Data Descriptions ***********************************************
        DataFile     % Name of the data file
        
        DataFreq     % Frequency of the data
               
        Data         % Data used in estimation

        FirstDay     % First day of the data
        
        LastDay      % Last day of the data
        
        T2M          % Time-to-maturities of the data
        
        % Estimation Spec *************************************************
        isCalib      % If 1, the model will be calibrated using
                     %     the parameter values stored in "ParamFile."
                     % Otherwise, the model will be estimated.
                         
        ParamsFile   % Name of the paramter file
        
        T2M2Use      % Maturities to be used in estimation
        
        nIterEKF     % If 0, the ordinary extended Kalman filter will be
                     %     implemented.
                     % If positive, the value represents # of iterations
                     %     in the iterated extende Kalman filter.
                     % If negative, the absolute value represents the torelance
                     %     size in the iterated extended Kalman filter.
                         
        dTau         % grid size when approximating forward rates
        
        MaxIter      % Max # of iterations in MATLAB optimazation routine
        
        isComputeSE  % If 1, std errors of the estimates will be computed.
        
        minDiff      % Min step size when computing Hessian matrix.
        
        EstStats     % Stats of fminisearch
        
        Time2Est     % estimation
        
        LLH          % log likelihood
    end
    
    %% B. Methods
    methods
        readData(obj) % Loads the input data
        
        convertDataFreq(obj) % Converts the data frequency
        
        setProperties(obj) % Loads and sets properties

        runEstimation(obj) % Runs estimation (filtering) routine
        
        [minusLLH, factor_mat] = runMLE(obj);  % Runs optimization
        
        [minusLLH, factor_mat] = runIEKF(obj) % Runs iterated extended Kalman filter

        [FR, varargout] = calcForwardRates(obj, factor, ELB, T2M) % Calcurates forward rates
        
        [R, varargout]  = calcSpotRates(obj, factor, ELB, T2M) % Calcurates spot rates

        Hessian      = calcHessian(obj) % Calcurates Hessian matrix
        
        yield_decomposition(obj, varargin) % decomposes yield changes
        
        WriteResult2Excel(obj, varargin) % saves the estimation results as
                                          % an Excel spreadsheet 
    end

end