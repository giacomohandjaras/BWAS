function replication_rate = compute_replication(bootstrap_output_discovery,bootstrap_output_replication,id_max_abs_corr)
% coded by Luca Cecchetti & Giacomo Handjaras, IMT-Lucca, Italy
% vers 20220619

% Timestamp
t0 = tic;

% Begin of the analysis part
fprintf('\n-----------------------------\n')
fprintf('| Estimate Replication Rate |\n')
fprintf('-----------------------------\n')

% Number of explored sample sizes
n_sample_sizes = size(bootstrap_output_discovery.simulated_corr,2);

% Number of actual correlations
n_actual_correlations = size(bootstrap_output_discovery.simulated_corr,3);

% Preallocate a matrix to store percentage of replication for uncorrected
% significant correlations
replication_percentage_uncorrected = ...
    nan(n_sample_sizes,n_actual_correlations);

% Preallocate a matrix to store percentage of replication for FDR corrected
% significant correlations
replication_percentage_fdr = ...
    nan(n_sample_sizes,n_actual_correlations);

% Preallocate a matrix to store percentage of replication for Bonferroni
% corrected significant correlations
replication_percentage_bonferroni = ...
    nan(n_sample_sizes,n_actual_correlations);

% For each sample size
for s = 1:n_sample_sizes
    
    % For each correlation
    for c = 1:n_actual_correlations
        
        % Find relationships passing the uncorrected significance threshold 
        % in the discovery sample
        temp_discovery_uncorrected = ...
            bootstrap_output_discovery.simulated_significant_pvalue_uncorrected(:,s,c);
        
        % Find relationships passing the uncorrected significance threshold 
        % in the replication sample        
        temp_replication_uncorrected = ...
            bootstrap_output_replication.simulated_significant_pvalue_uncorrected(:,s,c);
        
        % The percentage of successful replications is the number of times 
        % the correlation passed the same significance threshold in both 
        % the discovery and replication samples across resamplings. This 
        % value is then divided by the number of significant correlations
        % in the discovery dataset across bootstraps. Please note that in 
        % the output matrix, NaN means that a correlation was not
        % significant both in the discovery and in the replication (i.e.,
        % 0/0 = NaN).        
        replication_percentage_uncorrected(s,c) = ...
            sum(...
            all(...
            [temp_discovery_uncorrected, temp_replication_uncorrected],...
            2)) ...
            ./sum(temp_discovery_uncorrected).*100;
        
        % Find relationships passing the FDR corrected significance 
        % threshold in the discovery sample
        temp_discovery_fdr = ...
            bootstrap_output_discovery.simulated_significant_pvalue_fdr(:,s,c);
        
        % Find relationships passing the FDR corrected significance 
        % threshold in the replication sample
        temp_replication_fdr = ...
            bootstrap_output_replication.simulated_significant_pvalue_fdr(:,s,c);
        
        % Percentage of replication for FDR corrected correlations
        replication_percentage_fdr(s,c) = ...
            sum(...
            all(...
            [temp_discovery_fdr, temp_replication_fdr],...
            2)) ...
            ./sum(temp_discovery_fdr).*100;
        
        % Find relationships passing the Bonferroni corrected significance 
        % threshold in the discovery sample
        temp_discovery_bonferroni = ...
            bootstrap_output_discovery.simulated_significant_pvalue_bonferroni(:,s,c);
        
        % Find relationships passing the Bonferroni corrected significance 
        % threshold in the replication sample
        temp_replication_bonferroni = ...
            bootstrap_output_replication.simulated_significant_pvalue_bonferroni(:,s,c);

        % Percentage of replication for Bonferroni corrected correlations
        replication_percentage_bonferroni(s,c) = ...
            sum(...
            all(...
            [temp_discovery_bonferroni, temp_replication_bonferroni],...
            2)) ...
            ./sum(temp_discovery_bonferroni).*100;        
        
    end
    
end

% Replication percentage for the best uncorrected significant correlations
best_replication_percentage_uncorrected = ...
    replication_percentage_uncorrected(:,id_max_abs_corr);

% Replication percentage for the best FDR corrected significant 
% correlations
best_replication_percentage_fdr = ...
    replication_percentage_fdr(:,id_max_abs_corr);

% Replication percentage for the best Bonferroni corrected significant 
% correlations
best_replication_percentage_bonferroni = ...
    replication_percentage_bonferroni(:,id_max_abs_corr);

% Average replication percentage for uncorrected significant correlations
average_replication_percentage_uncorrected = ...
    nanmean(replication_percentage_uncorrected,2);

% Average replication percentage for FDR corrected correlations
average_replication_percentage_fdr = ...
    nanmean(replication_percentage_fdr,2);

% Average replication percentage for Bonferroni corrected correlations
average_replication_percentage_bonferroni = ...
    nanmean(replication_percentage_bonferroni,2);


% Prepare output

% Store replication rate for uncorrected significant correlations
replication_rate.uncorrected = replication_percentage_uncorrected;

% Store average replication rate for uncorrected significant correlations
replication_rate.average_uncorrected = ...
    average_replication_percentage_uncorrected;

% Store replication rate for the best uncorrected significant correlation
replication_rate.best_uncorrected = ...
    best_replication_percentage_uncorrected;

% Store replication rate for FDR corrected significant correlations
replication_rate.fdr = replication_percentage_fdr;

% Store average replication rate for FDR corrected significant correlations
replication_rate.average_fdr = ...
    average_replication_percentage_fdr;

% Store replication rate for the best FDR corrected significant correlation
replication_rate.best_fdr = ...
    best_replication_percentage_fdr;

% Store replication rate for Bonferroni significant correlations
replication_rate.bonferroni = replication_percentage_bonferroni;

% Store average replication rate for Bonferroni corrected significant 
% correlations
replication_rate.average_bonferroni = ...
    average_replication_percentage_bonferroni;

% Store replication rate for the best Bonferroni corrected significant 
% correlation
replication_rate.best_bonferroni = ...
    best_replication_percentage_bonferroni;

% Timestamp
elapsed = seconds(toc(t0));
elapsed.Format = 'hh:mm:ss.SSS';
elapsed = string(elapsed);

% Provide feedback
fprintf('All done! Elapsed time: %s\n',...
    elapsed);

end
