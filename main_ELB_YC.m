%% main_ELB_YC.m
%
%........................................................................
% Create: March 3, 2020 (Masaki Tanaka, Bank of Japan) 
%

%% **************************************************************
%% 0. Preamble
close all;
clear;
warning off MATLAB:xlswrite:AddSheet 
tic;

addpath('mttools');

%% **************************************************************
%% 1. Global options

nIterEst  = 1;     % # of iterations of estimation
Class2Use = 'BAFNS';
isEst     = 1;
File2Load = 'BAFNS_all_estimation_results.mat';

%% **************************************************************
%% 2. Model estimation
if isEst
    wtbar   = waitbar(0, 'Starting the estimation routine...', ...
        'Name', 'Runing MLE repeatedly... ');

    LLH_vec = zeros(nIterEst, 1);

    for i = 1 : nIterEst
        wtbar   = waitbar(i/(nIterEst + 1), wtbar,...
            sprintf( ['Implementing the MEL routine, ',...
            'changing the initial values (%d/%d) ...'], i, nIterEst));

        % Create a model instance
        eval(sprintf('instance_i = %s;', Class2Use));
    
        % Set Properties
        instance_i.setProperties;
    
        % Load the zero-coupon nominal yields
        instance_i.readData;
           
        % Load the parameter values
        params = readtable(instance_i.ParamsFile);
        if i > 1
            ELBplusObsErr = length(instance_i.T2M2Use) + instance_i.isELBmodel;
            chg_init_val = 1 + (rand(1, 19) - 0.5)/10;
            chg_init_val = [chg_init_val, ones(1, ELBplusObsErr)]; 
            params{1, :} = params{1, :} .* chg_init_val;
        end
        instance_i.Parameters = params;
    
        % Run the estimation routine    
        instance_i.runEstimation;
        LLH_vec(i) = instance_i.LLH;
    
        % Give the instance a unique name
        name_i = sprintf('Instance%d', i);
        eval(sprintf('%s = instance_i;', name_i));
    end

    [maxLLH, argmaxLLH] = max(LLH_vec);

    % Record the log
    LogName = sprintf('%s_%s_log.txt', datestr(datetime, 'yymmdd'), Class2Use);
    LogFile = fopen(LogName,'w');
    fprintf(LogFile, '-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-\n');
    fprintf(LogFile, '  Likelihood reaches at its maximum value in the trial %d.\n', argmaxLLH);
    fprintf(LogFile, '  The maximum value is %.3f.                              \n', maxLLH);
    fprintf(LogFile, '-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-\n');
    fclose(LogFile);

    % Save the estimation result
    eval(sprintf('save %s_all_estimation_results %s LLH_vec;',...
        Class2Use, sprintf('Instance%d ', (1:1:nIterEst)) )); 

else
    wtbar   = waitbar(1/2, 'Loading the existing estimation result...', ...
        'Name', 'Analyzing the estimated model... '); 
    load(File2Load);
    [maxLLH, argmaxLLH] = max(LLH_vec);
end
    
%% **************************************************************
%% 3. Analyze the estimated model
wtbar   = waitbar( 1 , wtbar, sprintf( ['Analyzing the estimated model ',...
        '(the result of trial No. %d) ...'], argmaxLLH));

File2Out = sprintf('%s_%s_result', datestr(datetime, 'yymmdd'), Class2Use);

% Copy and rename the model which has the the highest likelihood 
BestInstance = eval(sprintf('Instance%d', argmaxLLH));

% Decompose the changes in the yield curve on the week when the ELB shifts
BestInstance.yield_decomposition('T2M2Decomp', [0.5, 1:1:20]);

% Write out the result as an Excel file
BestInstance.WriteResult2Excel([File2Out, '.xlsx']);

% Save the analysis result as a matlab data file
save(File2Out, 'BestInstance'); 

%% **************************************************************
%% Postamble
fprintf('******************************************\n')
fprintf('   %s has successfuly finished!  \n', mfilename);
fprintf('******************************************\n')

record = toc;
close(wtbar);
delete(wtbar);
warning on MATLAB:xlswrite:AddSheet 