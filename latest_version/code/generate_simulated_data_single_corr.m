function simulated_data = generate_simulated_data_single_corr(corr_values, n_full_sample, tolerance)
% corr_value: correlation value in the full sample
%
% n_full_sample: the number of participants included in the full sample
% from which the reference ROIs are selected. As reported in the article,
% n_full_sample = 3928 for cortical thickness and n = 3604 for RSFC
%
% tolerance: a parameter to establish how much brain-behavior correlations
% derived from simulated data should differ from actual correlations. We
% set tolerance to 1e-4. Using smaller values could result in longer
% computation times
%
% Usage:
% [behavior, brain] = generate_simulated_data(0.2, ...
%     3604, ...
%     1e-4);
% coded by Luca Cecchetti & Giacomo Handjaras, IMT-Lucca, Italy
% vers 20220619

% Timestamp
t0 = tic;

% Begin of the analysis part
fprintf('\n---------------------------\n')
fprintf('| Generate Simulated Data |\n')
fprintf('---------------------------\n')

% Number of correlations (i.e., number of edges/ROIs/networks/components
% included in the analysis)
n_actual_correlations = numel(corr_values);

% Noise level. Please note that this parameter does not affect the
% correlation between brain and behavior. It is simply used to add noise to
% each of the two correlated measures, that's we did not provide an option
% to change it
sigma = .1;

% Generate random behavioral scores
behavior = sigma .* randn(n_full_sample, 1);

% Preallocate a matrix to store brain data correlated with behavior
brain = nan(n_full_sample, n_actual_correlations);

% An array to store correlation values for simulated data
simulated_correlations = nan(n_actual_correlations,1);

% For each correlation
for c = 1:n_actual_correlations
    
    % Initialize absolute rho difference to a large value so that for each
    % correlation we enter the while loop
    rho_diff = 1;
    
    % Until absolute rho difference between actual and simulated
    % correlation is larger than tolerance
    while rho_diff > tolerance
        
        % Generate random brain data
        brain_temp = sigma .* randn(n_full_sample, 1);
        
        % Extract the correlation coefficient of real data
        rho = corr_values(c);
        
        % Impose correlation between behavioral scores and brain data
        brain_temp_corr = ...
            (rho/abs(rho)) * ...
            sqrt(rho^2) * ...
            behavior + ...
            sqrt(1-rho^2) * ...
            brain_temp;
        
        % Get the correlation between simulated brain and behavior
        synthetic_rho = corr(behavior,brain_temp_corr);
        
        % Estimate the absolute difference between actual and simulated
        % correlations
        rho_diff = abs(rho-synthetic_rho);
        
    end
    
    % Store simulated brain data
    brain(:,c) = brain_temp_corr;
    
    % Store correlation for simulated data
    simulated_correlations(c) = synthetic_rho;    
    
end

% Prepare output structure

% Behavior in the full sample
simulated_data.behavior = behavior;

% Brain features in the full sample
simulated_data.brain = brain;

% Original correlation values
simulated_data.original_corr_values = corr_values;

% Correlation values for simulated data
simulated_data.simulated_corr_values = simulated_correlations;

% Timestamp
elapsed = seconds(toc(t0));
elapsed.Format = 'hh:mm:ss.SSS';
elapsed = string(elapsed);

% Provide feedback
fprintf('All done! Elapsed time: %s\n',...
    elapsed);


end
