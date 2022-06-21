function simulated_data = generate_simulated_data(source_data, n_full_sample, tolerance)
% source_data: the path to a csv file storing correlation values from the
% Marek, Tervo-Clemmens et al. 2022 article. Source data are publicly
% available and can be downloaded from the article webpage. The number of
% rows in the csv relates to the number of edges/ROIs/components/networks
% included in the analysis. For instance: components_cognition.csv is a
% text file storing correlation values between the 100 components derived
% from RSFC analysis and cognitive abilities (Fig. 1B left panel)
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
% [behavior, brain] = generate_simulated_data('components_cognition.csv', ...
%     3604, ...
%     1e-4);
% coded by Luca Cecchetti & Giacomo Handjaras, IMT-Lucca, Italy
% vers 20220619

% Timestamp
t0 = tic;

% Number of steps in progress bar
step_progress = 10;

% Begin of the analysis part
fprintf('\n---------------------------\n')
fprintf('| Generate Simulated Data |\n')
fprintf('---------------------------\n')

% Read correlations in csv file
actual_correlations = csvread(source_data);

% Number of correlations (i.e., number of edges/ROIs/networks/components
% included in the analysis)
n_actual_correlations = numel(actual_correlations);

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
        rho = actual_correlations(c);
        
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
    
    % Print feedback about progress
    if mod(c,floor(n_actual_correlations/step_progress))==0
        
        % Time stamp
        elapsed = seconds(toc(t0));
        elapsed.Format = 'hh:mm:ss.SSS';
        elapsed = string(elapsed);
        
        % Provide feedback
        fprintf('Simulated %d out of %d correlations - Elapsed time: %s\n',...
            c,n_actual_correlations,elapsed);
        
    end
    
    
end

% Number of observations in the discovery dataset
n_observations_discovery = round(n_full_sample/2);

% Generate an array of values to randomly split the original full sample in
% discovery and replication
split_half_random = randperm(n_full_sample);

% Randomly attribute behavioral scores of half of the participants to the
% discovery dataset
behavior_discovery = ...
    behavior(split_half_random(1:n_observations_discovery));

% Attribute the other half of the participants behavioral scores to the
% replication dataset
behavior_replication = ...
    behavior(split_half_random(n_observations_discovery+1:end));

% Randomly attribute brain data of half of the participants to the
% discovery dataset
brain_discovery = ...
    brain(split_half_random(1:n_observations_discovery),:);

% Attribute the other half of the participants brain data to the
% replication dataset
brain_replication = ...
    brain(split_half_random(n_observations_discovery+1:end),:);

% Prepare output structure

% Behavior in the full sample
simulated_data.behavior = behavior;

% Brain features in the full sample
simulated_data.brain = brain;

% Behavior in the discovery sample
simulated_data.discovery_behavior = behavior_discovery;

% Behavior in the replication sample
simulated_data.replication_behavior = behavior_replication;

% Brain in the discovery sample
simulated_data.discovery_brain = brain_discovery;

% Brain in the replication sample
simulated_data.replication_brain = brain_replication;

% Correlation values in the original Figure 1
simulated_data.original_corr_values = actual_correlations;

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
