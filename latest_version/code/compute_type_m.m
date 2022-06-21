function type_m_errors = compute_type_m(bootstrap_output)
% coded by Luca Cecchetti & Giacomo Handjaras, IMT-Lucca, Italy
% vers 20220619

% Timestamp
t0 = tic;

% Begin of the analysis part
fprintf('\n-------------------------\n')
fprintf('| Estimate Type M error |\n')
fprintf('-------------------------\n')

% Number of explored sample sizes
n_sample_sizes = size(bootstrap_output.simulated_corr,2);

% Uncorrected significance at each sample size and resampling for
% correlations passing threshold for uncorrected significance in the
% full sample
significant_uncorrected = ...
    bootstrap_output.simulated_significant_pvalue_uncorrected(...
    :,:,bootstrap_output.p_full_significant_uncorrected);

% FDR corrected significance at each sample size and resampling for
% correlations passing threshold for FDR corrected significance in the
% full sample
significant_fdr = ...
    bootstrap_output.simulated_significant_pvalue_fdr(...
    :,:,bootstrap_output.p_full_significant_fdr);

% Bonferroni corrected significance at each sample size and resampling for
% correlations passing threshold for Bonferroni corrected significance in
% the full sample
significant_bonferroni = ...
    bootstrap_output.simulated_significant_pvalue_bonferroni(...
    :,:,bootstrap_output.p_full_significant_bonferroni);

% Uncorrected significant correlations at each sample size and resampling
% for correlations passing threshold for uncorrected significance in the
% full sample
significant_uncorrected_corr = ...
    bootstrap_output.simulated_corr(...
    :,:,bootstrap_output.p_full_significant_uncorrected);

% FDR corrected significant correlations at each sample size and resampling
% for correlations passing threshold for FDR corrected significance in the
% full sample
significant_fdr_corr = ...
    bootstrap_output.simulated_corr(...
    :,:,bootstrap_output.p_full_significant_fdr);

% Bonferroni corrected significant correlations at each sample size and
% resampling for correlations passing threshold for FDR corrected
% significance in the full sample
significant_bonferroni_corr = ...
    bootstrap_output.simulated_corr(...
    :,:,bootstrap_output.p_full_significant_bonferroni);

% Uncorrected significant correlations in the full sample
significant_uncorrected_corr_full = ...
    bootstrap_output.r_full(...
    bootstrap_output.p_full_significant_uncorrected);

% FDR corrected significant correlations in the full sample
significant_fdr_corr_full = ...
    bootstrap_output.r_full(...
    bootstrap_output.p_full_significant_fdr);

% Bonferroni corrected significant correlations in the full sample
significant_bonferroni_corr_full = ...
    bootstrap_output.r_full(...
    bootstrap_output.p_full_significant_bonferroni);

% Number of correlations passing the uncorrected threshold for significance
% in the full sample
n_actual_correlations_uncorrected = size(significant_uncorrected,3);

% Number of correlations passing the FDR corrected threshold for
% significance in the full sample
n_actual_correlations_fdr = size(significant_fdr,3);

% Number of correlations passing the Bonferroni corrected threshold for
% significance in the full sample
n_actual_correlations_bonferroni = size(significant_bonferroni,3);

% Number of correlations passing the uncorrected significance threshold and
% inflated by 50% with respect to the full sample
inflation_50_uncorrected = ...
    nan(n_sample_sizes,n_actual_correlations_uncorrected);

% Number of correlations passing the FDR corrected significance threshold
% and inflated by 50% with respect to the full sample
inflation_50_fdr = ...
    nan(n_sample_sizes,n_actual_correlations_fdr);

% Number of correlations passing the Bonferroni corrected significance
% threshold and inflated by 50% with respect to the full sample
inflation_50_bonferroni = ...
    nan(n_sample_sizes,n_actual_correlations_bonferroni);

% Number of correlations passing the uncorrected significance threshold and
% inflated by 100% with respect to the full sample
inflation_100_uncorrected = ...
    nan(n_sample_sizes,n_actual_correlations_uncorrected);

% Number of correlations passing the FDR corrected significance threshold
% and inflated by 100% with respect to the full sample
inflation_100_fdr = ...
    nan(n_sample_sizes,n_actual_correlations_fdr);

% Number of correlations passing the Bonferroni corrected significance
% threshold and inflated by 100% with respect to the full sample
inflation_100_bonferroni = ...
    nan(n_sample_sizes,n_actual_correlations_bonferroni);

% Number of correlations passing the uncorrected significance threshold and
% inflated by 200% with respect to the full sample
inflation_200_uncorrected = ...
    nan(n_sample_sizes,n_actual_correlations_uncorrected);

% Number of correlations passing the FDR corrected significance threshold
% and inflated by 200% with respect to the full sample
inflation_200_fdr = ...
    nan(n_sample_sizes,n_actual_correlations_fdr);

% Number of correlations passing the Bonferroni corrected significance
% threshold  and inflated by 200% with respect to the full sample
inflation_200_bonferroni = ...
    nan(n_sample_sizes,n_actual_correlations_bonferroni);

% For each sample size
for s = 1:n_sample_sizes
    
    % For each correlation passing the uncorrected significance threshold
    % in the full sample
    for u = 1:n_actual_correlations_uncorrected
        
        % Get correlations of all resamplings
        temp_correlations = ...
            significant_uncorrected_corr(...
            significant_uncorrected(:,s,u),s,u);
        
        % Compute inflation rate dividing the correlations found across
        % resamplings by the correlation in the full sample. Each element
        % in the inflation rate array expresses the magnitude of inflation:
        % >1.5 is >50% inflation, >2 is >100% inflation, and >3 is >200%
        % inflation.
        inflation_rate = ...
            temp_correlations./significant_uncorrected_corr_full(u);
        
        % Remove negative values because they indicate that correlation in
        % one resampling has opposite sign with respect to the full sample
        inflation_rate = inflation_rate(inflation_rate>0);
        
        % Percentage of correlations inflated by 50%
        inflation_50_uncorrected(s,u) = ...
            sum(inflation_rate>1.5)/numel(inflation_rate)*100;
        
        % Percentage of correlations inflated by 100%
        inflation_100_uncorrected(s,u) = ...
            sum(inflation_rate>2)/numel(inflation_rate)*100;
        
        % Percentage of correlations inflated by 200%
        inflation_200_uncorrected(s,u) = ...
            sum(inflation_rate>3)/numel(inflation_rate)*100;
        
    end
    
    % For each correlation passing the FDR corrected significance threshold
    % in the full sample
    for f = 1:n_actual_correlations_fdr
        
        % Get correlations of all resamplings
        temp_correlations = ...
            significant_fdr_corr(...
            significant_fdr(:,s,f),s,f);
        
        % Compute inflation rate dividing the correlations found across
        % resamplings by the correlation in the full sample. Each element
        % in the inflation rate array expresses the magnitude of inflation:
        % >1.5 is >50% inflation, >2 is >100% inflation, and >3 is >200%
        % inflation.
        inflation_rate = ...
            temp_correlations./significant_fdr_corr_full(f);
        
        % Remove negative values because they indicate that correlation in
        % one resampling has opposite sign with respect to the full sample
        inflation_rate = inflation_rate(inflation_rate>0);
        
        % Percentage of correlations inflated by 50%
        inflation_50_fdr(s,f) = ...
            sum(inflation_rate>1.5)/numel(inflation_rate)*100;
        
        % Percentage of correlations inflated by 100%
        inflation_100_fdr(s,f) = ...
            sum(inflation_rate>2)/numel(inflation_rate)*100;
        
        % Percentage of correlations inflated by 200%
        inflation_200_fdr(s,f) = ...
            sum(inflation_rate>3)/numel(inflation_rate)*100;
        
    end
    
    % For each correlation passing the Bonferroni corrected significance
    % threshold in the full sample
    for b = 1:n_actual_correlations_bonferroni
        
        % Get correlations of all resamplings
        temp_correlations = ...
            significant_bonferroni_corr(...
            significant_bonferroni(:,s,b),s,b);
        
        % Compute inflation rate dividing the correlations found across
        % resamplings by the correlation in the full sample. Each element
        % in the inflation rate array expresses the magnitude of inflation:
        % >1.5 is >50% inflation, >2 is >100% inflation, and >3 is >200%
        % inflation.
        inflation_rate = ...
            temp_correlations./significant_bonferroni_corr_full(b);
        
        % Remove negative values because they indicate that correlation in
        % one resampling has opposite sign with respect to the full sample
        inflation_rate = inflation_rate(inflation_rate>0);
        
        % Percentage of correlations inflated by 50%
        inflation_50_bonferroni(s,b) = ...
            sum(inflation_rate>1.5)/numel(inflation_rate)*100;
        
        % Percentage of correlations inflated by 100%
        inflation_100_bonferroni(s,b) = ...
            sum(inflation_rate>2)/numel(inflation_rate)*100;
        
        % Percentage of correlations inflated by 200%
        inflation_200_bonferroni(s,b) = ...
            sum(inflation_rate>3)/numel(inflation_rate)*100;
        
    end
    
end

% Compute in each sample size the percentage of inflation by 50% for 
% uncorrected significance averaged across correlations
average_inflation_50_uncorrected = ...
    nanmean(inflation_50_uncorrected,2);

% Compute in each sample size the percentage of inflation by 100% for 
% uncorrected significance averaged across correlations
average_inflation_100_uncorrected = ...
    nanmean(inflation_100_uncorrected,2);

% Compute in each sample size the percentage of inflation by 200% for 
% uncorrected significance averaged across correlations
average_inflation_200_uncorrected = ...
    nanmean(inflation_200_uncorrected,2);

% Compute in each sample size the percentage of inflation by 50% for 
% FDR corrected significance averaged across correlations
average_inflation_50_fdr = ...
    nanmean(inflation_50_fdr,2);

% Compute in each sample size the percentage of inflation by 100% for 
% FDR corrected significance averaged across correlations
average_inflation_100_fdr = ...
    nanmean(inflation_100_fdr,2);

% Compute in each sample size the percentage of inflation by 200% for 
% FDR corrected significance averaged across correlations
average_inflation_200_fdr = ...
    nanmean(inflation_200_fdr,2);

% Compute in each sample size the percentage of inflation by 50% for 
% Bonferroni corrected significance averaged across correlations
average_inflation_50_bonferroni = ...
    nanmean(inflation_50_bonferroni,2);

% Compute in each sample size the percentage of inflation by 100% for 
% Bonferroni corrected significance averaged across correlations
average_inflation_100_bonferroni = ...
    nanmean(inflation_100_bonferroni,2);

% Compute in each sample size the percentage of inflation by 200% for 
% Bonferroni corrected significance averaged across correlations
average_inflation_200_bonferroni = ...
    nanmean(inflation_200_bonferroni,2);

% Prepare Output

% Store 50% inflation for the uncorrected threshold
type_m_errors.inflation_50_uncorrected = inflation_50_uncorrected;

% Store average 50% inflation for the uncorrected threshold
type_m_errors.average_inflation_50_uncorrected = ...
    average_inflation_50_uncorrected;

% Store 100% inflation for the uncorrected threshold
type_m_errors.inflation_100_uncorrected = inflation_100_uncorrected;

% Store average 100% inflation for the uncorrected threshold
type_m_errors.average_inflation_100_uncorrected = ...
    average_inflation_100_uncorrected;

% Store 200% inflation for the uncorrected threshold
type_m_errors.inflation_200_uncorrected = inflation_200_uncorrected;

% Store average 200% inflation for the uncorrected threshold
type_m_errors.average_inflation_200_uncorrected = ...
    average_inflation_200_uncorrected;

% Store 50% inflation for the FDR corrected threshold
type_m_errors.inflation_50_fdr = inflation_50_fdr;

% Store average 50% inflation for the FDR corrected threshold
type_m_errors.average_inflation_50_fdr = ...
    average_inflation_50_fdr;

% Store 100% inflation for the FDR corrected threshold
type_m_errors.inflation_100_fdr = inflation_100_fdr;

% Store average 100% inflation for the FDR corrected threshold
type_m_errors.average_inflation_100_fdr = ...
    average_inflation_100_fdr;

% Store 200% inflation for the FDR corrected threshold
type_m_errors.inflation_200_fdr = inflation_200_fdr;

% Store average 200% inflation for the FDR corrected threshold
type_m_errors.average_inflation_200_fdr = ...
    average_inflation_200_fdr;

% Store 50% inflation for the Bonferroni corrected threshold
type_m_errors.inflation_50_bonferroni = inflation_50_bonferroni;

% Store average 50% inflation for the Bonferroni corrected threshold
type_m_errors.average_inflation_50_bonferroni = ...
    average_inflation_50_bonferroni;

% Store 100% inflation for the Bonferroni corrected threshold
type_m_errors.inflation_100_bonferroni = inflation_100_bonferroni;

% Store average 100% inflation for the Bonferroni corrected threshold
type_m_errors.average_inflation_100_bonferroni = ...
    average_inflation_100_bonferroni;

% Store 200% inflation for the Bonferroni corrected threshold
type_m_errors.inflation_200_bonferroni = inflation_200_bonferroni;

% Store average 200% inflation for the Bonferroni corrected threshold
type_m_errors.average_inflation_200_bonferroni = ...
    average_inflation_200_bonferroni;

% Timestamp
elapsed = seconds(toc(t0));
elapsed.Format = 'hh:mm:ss.SSS';
elapsed = string(elapsed);

% Provide feedback
fprintf('All done! Elapsed time: %s\n',...
    elapsed);

end





