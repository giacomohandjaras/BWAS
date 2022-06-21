function type_m_errors = compute_type_m_alternative(bootstrap_output,aggregation_method)
% coded by Luca Cecchetti & Giacomo Handjaras, IMT-Lucca, Italy
% vers 20220619

% To assess type M error, we computed the ratio between the correlation value obtained 
% for each significant resampling and the one obtained in the full sample. 
% The inflation percentage was then estimated at each sample size as the average across resamplings.
% This procedure is different from the one proposed in M&TC 
% (use our compute_type_m.m to make a direct comparison with  M&TC).
% M&TC determined the percentage of results that were inflated, 
% using multiple bins of varying magnitudes of inflation (50%, 100% and 200%).

% Timestamp
t0 = tic;

% Begin of the analysis part
fprintf('\n-------------------------\n')
fprintf('| Estimate Type M error |\n')
fprintf('-------------------------\n')

% The best correlation in the full sample
[~, id_max_abs_corr] = max(abs(bootstrap_output.r_full));

% Position of the best correlation among those passing the uncorrected
% significance threshold
id_max_abs_corr_uncorrected = ...
    sum(bootstrap_output.p_full_significant_uncorrected(...
    1:id_max_abs_corr));

% Position of the best correlation among those passing the FDR corrected
% significance threshold
id_max_abs_corr_fdr = ...
    sum(bootstrap_output.p_full_significant_fdr(...
    1:id_max_abs_corr));

% Position of the best correlation among those passing the Bonferroni
% corrected significance threshold
id_max_abs_corr_bonferroni = ...
    sum(bootstrap_output.p_full_significant_bonferroni(...
    1:id_max_abs_corr));

% Number of explored sample sizes
n_sample_sizes = size(bootstrap_output.simulated_corr,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UNCORRECTED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Uncorrected significance at each sample size and resampling for
% correlations passing threshold for uncorrected significance in the
% full sample
significant_uncorrected = ...
    bootstrap_output.simulated_significant_pvalue_uncorrected(...
    :,:,bootstrap_output.p_full_significant_uncorrected);

% Correlations at each sample size and resampling passing threshold for
% uncorrected significance in the full sample
significant_uncorrected_corr = ...
    bootstrap_output.simulated_corr(...
    :,:,bootstrap_output.p_full_significant_uncorrected);

% Uncorrected significant correlations in the full sample
significant_uncorrected_corr_full = ...
    bootstrap_output.r_full(...
    bootstrap_output.p_full_significant_uncorrected);

% Number of correlations passing the uncorrected threshold for significance
% in the full sample
n_actual_correlations_uncorrected = size(significant_uncorrected,3);

% Preallocate an array to store average inflation percentage for each
% correlation passing the uncorrected threshold at all sample sizes
aggregated_inflation_percentage_thresholded_uncorrected = ...
    nan(n_sample_sizes,n_actual_correlations_uncorrected);

% Preallocate an array to store average inflation percentage for each
% correlation passing the uncorrected threshold in the full sample,
% regardless of significance at each resampling
aggregated_inflation_percentage_uncorrected = ...
    nan(n_sample_sizes,n_actual_correlations_uncorrected);

% For each sample size
for s = 1:n_sample_sizes
    
    % For each correlation passing the uncorrected significance threshold
    % in the full sample
    for u = 1:n_actual_correlations_uncorrected
        
        % Get correlations of all resamplings passing the significance both
        % in the full sample and resampling
        temp_correlations_thresholded = ...
            significant_uncorrected_corr(...
            significant_uncorrected(:,s,u),s,u);
        
        % Get correlations of all resamplings passing the significance in
        % the full sample
        temp_correlations = significant_uncorrected_corr(:,s,u);
        
        % Compute inflation rate dividing the correlations found across
        % resamplings by the correlation in the full sample. Each element
        % in the inflation rate array expresses the magnitude of inflation:
        % 1.5 is 50% inflation, 2 is 100% inflation, and 3 is 200%
        % inflation.
        inflation_rate_thresholded = ...
            temp_correlations_thresholded ./ ...
            significant_uncorrected_corr_full(u);
        
        % Compute inflation rate dividing the correlations found across
        % resamplings by the correlation in the full sample. Each element
        % in the inflation rate array expresses the magnitude of inflation:
        % 1.5 is 50% inflation, 2 is 100% inflation, and 3 is 200%
        % inflation.
        inflation_rate = ...
            temp_correlations ./ ...
            significant_uncorrected_corr_full(u);
        
        % Remove negative values because they indicate that correlation in
        % one resampling has opposite sign with respect to the full sample
        inflation_rate_thresholded = ...
            inflation_rate_thresholded(inflation_rate_thresholded > 0);
        
        % Remove negative values because they indicate that correlation in
        % one resampling has opposite sign with respect to the full sample
        inflation_rate = ...
            inflation_rate(inflation_rate > 0);
        
        % Transform inflation in percentage. Negative values mean shrinking
        % with respect to the correlation in the full sample
        inflation_percentage_thresholded = ...
            (inflation_rate_thresholded-1).*100;
        
        % Transform inflation in percentage. Negative values mean shrinking
        % with respect to the correlation in the full sample
        inflation_percentage = ...
            (inflation_rate-1).*100;
        
        % Aggregated inflation across resamplings
        switch lower(aggregation_method)
            
            case 'mean'
                
                % Compute the mean inflation percentage across resamplings
                % for correlations passing the uncorrected threshold in
                % both full sample and resampling
                aggregated_inflation_percentage_thresholded_uncorrected(s,u) = ...
                    mean(inflation_percentage_thresholded);
                
                % Compute the mean inflation percentage across resamplings
                % for correlations passing the uncorrected threshold in
                % the full sample
                aggregated_inflation_percentage_uncorrected(s,u) = ...
                    mean(inflation_percentage);
                
            case 'median'
                
                % Compute the median inflation percentage across
                % resamplings for correlations passing the uncorrected
                % threshold in both full sample and resampling
                aggregated_inflation_percentage_thresholded_uncorrected(s,u) = ...
                    median(inflation_percentage_thresholded);
                
                % Compute the median inflation percentage across
                % resamplings for correlations passing the uncorrected
                % threshold in the full sample
                aggregated_inflation_percentage_uncorrected(s,u) = ...
                    median(inflation_percentage);
                
        end
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FDR CORRECTED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FDR corrected significance at each sample size and resampling for
% correlations passing threshold for FDR corrected significance in the
% full sample
significant_fdr = ...
    bootstrap_output.simulated_significant_pvalue_fdr(...
    :,:,bootstrap_output.p_full_significant_fdr);

% Correlations at each sample size and resampling passing threshold for
% FDR corrected significance in the full sample
significant_fdr_corr = ...
    bootstrap_output.simulated_corr(...
    :,:,bootstrap_output.p_full_significant_fdr);

% FDR corrected significant correlations in the full sample
significant_fdr_corr_full = ...
    bootstrap_output.r_full(...
    bootstrap_output.p_full_significant_fdr);

% Number of correlations passing the FDR corrected threshold for
% significance in the full sample
n_actual_correlations_fdr = size(significant_fdr,3);

% Preallocate an array to store average inflation percentage for each
% correlation passing the FDR corrected threshold at all sample sizes
aggregated_inflation_percentage_thresholded_fdr = ...
    nan(n_sample_sizes,n_actual_correlations_fdr);

% Preallocate an array to store average inflation percentage for each
% correlation passing the FDR corrected threshold in the full sample,
% regardless of significance at each resampling
aggregated_inflation_percentage_fdr = ...
    nan(n_sample_sizes,n_actual_correlations_fdr);

% For each sample size
for s = 1:n_sample_sizes
    
    % For each correlation passing the FDR corrected significance threshold
    % in the full sample
    for u = 1:n_actual_correlations_fdr
        
        % Get correlations of all resamplings passing the significance both
        % in the full sample and resampling
        temp_correlations_thresholded = ...
            significant_fdr_corr(...
            significant_fdr(:,s,u),s,u);
        
        % Get correlations of all resamplings passing the significance in
        % the full sample
        temp_correlations = significant_fdr_corr(:,s,u);
        
        % Compute inflation rate dividing the correlations found across
        % resamplings by the correlation in the full sample. Each element
        % in the inflation rate array expresses the magnitude of inflation:
        % 1.5 is 50% inflation, 2 is 100% inflation, and 3 is 200%
        % inflation.
        inflation_rate_thresholded = ...
            temp_correlations_thresholded ./ ...
            significant_fdr_corr_full(u);
        
        % Compute inflation rate dividing the correlations found across
        % resamplings by the correlation in the full sample. Each element
        % in the inflation rate array expresses the magnitude of inflation:
        % 1.5 is 50% inflation, 2 is 100% inflation, and 3 is 200%
        % inflation.
        inflation_rate = ...
            temp_correlations ./ ...
            significant_fdr_corr_full(u);
        
        % Remove negative values because they indicate that correlation in
        % one resampling has opposite sign with respect to the full sample
        inflation_rate_thresholded = ...
            inflation_rate_thresholded(inflation_rate_thresholded > 0);
        
        % Remove negative values because they indicate that correlation in
        % one resampling has opposite sign with respect to the full sample
        inflation_rate = ...
            inflation_rate(inflation_rate > 0);
        
        % Transform inflation in percentage. Negative values mean shrinking
        % with respect to the correlation in the full sample
        inflation_percentage_thresholded = ...
            (inflation_rate_thresholded-1).*100;
        
        % Transform inflation in percentage. Negative values mean shrinking
        % with respect to the correlation in the full sample
        inflation_percentage = ...
            (inflation_rate-1).*100;
        
        % Aggregated inflation across resamplings
        switch lower(aggregation_method)
            
            case 'mean'
                
                % Compute the mean inflation percentage across resamplings
                % for correlations passing the FDR corrected threshold in
                % both full sample and resampling
                aggregated_inflation_percentage_thresholded_fdr(s,u) = ...
                    mean(inflation_percentage_thresholded);
                
                % Compute the mean inflation percentage across resamplings
                % for correlations passing the FDR corrected threshold in
                % the full sample
                aggregated_inflation_percentage_fdr(s,u) = ...
                    mean(inflation_percentage);
                
            case 'median'
                
                % Compute the median inflation percentage across
                % resamplings for correlations passing the FDR corrected
                % threshold in both full sample and resampling
                aggregated_inflation_percentage_thresholded_fdr(s,u) = ...
                    median(inflation_percentage_thresholded);
                
                % Compute the median inflation percentage across
                % resamplings for correlations passing the FDR corrected
                % threshold in the full sample
                aggregated_inflation_percentage_fdr(s,u) = ...
                    median(inflation_percentage);
                
        end
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% BONFERRONI CORRECTED %%%%%%%%%%%%%%%%%%%%%%%%%%

% Bonferroni corrected significance at each sample size and resampling for
% correlations passing threshold for Bonferroni corrected significance in
% the full sample
significant_bonferroni = ...
    bootstrap_output.simulated_significant_pvalue_bonferroni(...
    :,:,bootstrap_output.p_full_significant_bonferroni);

% Correlations at each sample size and resampling passing Bonferroni
% threshold for significance in the full sample
significant_bonferroni_corr = ...
    bootstrap_output.simulated_corr(...
    :,:,bootstrap_output.p_full_significant_bonferroni);

% Bonferroni corrected significant correlations in the full sample
significant_bonferroni_corr_full = ...
    bootstrap_output.r_full(...
    bootstrap_output.p_full_significant_bonferroni);

% Number of correlations passing the Bonferroni corrected threshold for
% significance in the full sample
n_actual_correlations_bonferroni = size(significant_bonferroni,3);

% Preallocate an array to store average inflation percentage for each
% correlation passing the Bonferroni corrected threshold at all sample
% sizes
aggregated_inflation_percentage_thresholded_bonferroni = ...
    nan(n_sample_sizes,n_actual_correlations_bonferroni);

% Preallocate an array to store average inflation percentage for each
% correlation passing the Bonferroni corrected threshold in the full
% sample, regardless of significance at each resampling
aggregated_inflation_percentage_bonferroni = ...
    nan(n_sample_sizes,n_actual_correlations_bonferroni);

% For each sample size
for s = 1:n_sample_sizes
    
    % For each correlation passing the Bonferroni corrected significance
    % threshold in the full sample
    for u = 1:n_actual_correlations_bonferroni
        
        % Get correlations of all resamplings passing the significance both
        % in the full sample and resampling
        temp_correlations_thresholded = ...
            significant_bonferroni_corr(...
            significant_bonferroni(:,s,u),s,u);
        
        % Get correlations of all resamplings passing the significance in
        % the full sample
        temp_correlations = significant_bonferroni_corr(:,s,u);
        
        % Compute inflation rate dividing the correlations found across
        % resamplings by the correlation in the full sample. Each element
        % in the inflation rate array expresses the magnitude of inflation:
        % 1.5 is 50% inflation, 2 is 100% inflation, and 3 is 200%
        % inflation.
        inflation_rate_thresholded = ...
            temp_correlations_thresholded ./ ...
            significant_bonferroni_corr_full(u);
        
        % Compute inflation rate dividing the correlations found across
        % resamplings by the correlation in the full sample. Each element
        % in the inflation rate array expresses the magnitude of inflation:
        % 1.5 is 50% inflation, 2 is 100% inflation, and 3 is 200%
        % inflation.
        inflation_rate = ...
            temp_correlations ./ ...
            significant_bonferroni_corr_full(u);
        
        % Remove negative values because they indicate that correlation in
        % one resampling has opposite sign with respect to the full sample
        inflation_rate_thresholded = ...
            inflation_rate_thresholded(inflation_rate_thresholded > 0);
        
        % Remove negative values because they indicate that correlation in
        % one resampling has opposite sign with respect to the full sample
        inflation_rate = ...
            inflation_rate(inflation_rate > 0);
        
        % Transform inflation in percentage. Negative values mean shrinking
        % with respect to the correlation in the full sample
        inflation_percentage_thresholded = ...
            (inflation_rate_thresholded-1).*100;
        
        % Transform inflation in percentage. Negative values mean shrinking
        % with respect to the correlation in the full sample
        inflation_percentage = ...
            (inflation_rate-1).*100;
        
        % Aggregated inflation across resamplings
        switch lower(aggregation_method)
            
            case 'mean'
                
                % Compute the mean inflation percentage across resamplings
                % for correlations passing the Bonferroni corrected
                % threshold in both full sample and resampling
                aggregated_inflation_percentage_thresholded_bonferroni(s,u) = ...
                    mean(inflation_percentage_thresholded);
                
                % Compute the mean inflation percentage across resamplings
                % for correlations passing the Bonferroni corrected
                % threshold in the full sample
                aggregated_inflation_percentage_bonferroni(s,u) = ...
                    mean(inflation_percentage);
                
            case 'median'
                
                % Compute the median inflation percentage across
                % resamplings for correlations passing the Bonferroni
                % corrected threshold in both full sample and resampling
                aggregated_inflation_percentage_thresholded_bonferroni(s,u) = ...
                    median(inflation_percentage_thresholded);
                
                % Compute the median inflation percentage across
                % resamplings for correlations passing the Bonferroni
                % corrected threshold in the full sample
                aggregated_inflation_percentage_bonferroni(s,u) = ...
                    median(inflation_percentage);
                
        end
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UNCORRECTED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mean inflation percentage across correlations passing the uncorrected
% significance threshold both in the full sample and in resamplings
average_inflation_thresholded_uncorrected = ...
    nanmean(aggregated_inflation_percentage_thresholded_uncorrected,2);

% Mean inflation percentage across correlations passing the uncorrected
% significance in the full sample
average_inflation_uncorrected = ...
    nanmean(aggregated_inflation_percentage_uncorrected,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FDR CORRECTED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mean inflation percentage across correlations passing the FDR corrected
% significance threshold both in the full sample and in resamplings
average_inflation_thresholded_fdr = ...
    nanmean(aggregated_inflation_percentage_thresholded_fdr,2);

% Mean inflation percentage across correlations passing the FDR corrected
% significance in the full sample
average_inflation_fdr = ...
    nanmean(aggregated_inflation_percentage_fdr,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%% BONFERRONI CORRECTED %%%%%%%%%%%%%%%%%%%%%%%%%%

% Mean inflation percentage across correlations passing the Bonferroni
% corrected significance threshold both in the full sample and in 
% resamplings
average_inflation_thresholded_bonferroni = ...
    nanmean(aggregated_inflation_percentage_thresholded_bonferroni,2);

% Mean inflation percentage across correlations passing the Bonferroni
% corrected significance in the full sample
average_inflation_bonferroni = ...
    nanmean(aggregated_inflation_percentage_bonferroni,2);

% Prepare output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UNCORRECTED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inflation percentage of thresholded uncorrected correlations
type_m_errors.inflation_percentage_thresholded_uncorrected = ...
    aggregated_inflation_percentage_thresholded_uncorrected;

% Inflation percentage of uncorrected correlations
type_m_errors.inflation_percentage_uncorrected = ...
    aggregated_inflation_percentage_uncorrected;

% Average inflation percentage of thresholded uncorrected correlations
type_m_errors.average_inflation_percentage_thresholded_uncorrected = ...
    average_inflation_thresholded_uncorrected;

% Average inflation percentage of uncorrected correlations
type_m_errors.average_inflation_percentage_uncorrected = ...
    average_inflation_uncorrected;

% Inflation percentage of the best correlation thresholded uncorrected
type_m_errors.best_inflation_percentage_thresholded_uncorrected = ...
    aggregated_inflation_percentage_thresholded_uncorrected(:,...
    id_max_abs_corr_uncorrected);

% Inflation percentage of the best correlation uncorrected
type_m_errors.best_inflation_percentage_uncorrected = ...
    aggregated_inflation_percentage_uncorrected(:,...
    id_max_abs_corr_uncorrected);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FDR CORRECTED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inflation percentage of thresholded FDR corrected correlations
type_m_errors.inflation_percentage_thresholded_fdr = ...
    aggregated_inflation_percentage_thresholded_fdr;

% Inflation percentage of FDR corrected correlations
type_m_errors.inflation_percentage_fdr = ...
    aggregated_inflation_percentage_fdr;

% Average inflation percentage of thresholded FDR corrected correlations
type_m_errors.average_inflation_percentage_thresholded_fdr = ...
    average_inflation_thresholded_fdr;

% Average inflation percentage of FDR corrected correlations
type_m_errors.average_inflation_percentage_fdr = ...
    average_inflation_fdr;

% Inflation percentage of the best correlation thresholded FDR corrected
type_m_errors.best_inflation_percentage_thresholded_fdr = ...
    aggregated_inflation_percentage_thresholded_fdr(:,...
    id_max_abs_corr_fdr);

% Inflation percentage of the best correlation FDR corrected
type_m_errors.best_inflation_percentage_fdr = ...
    aggregated_inflation_percentage_fdr(:,...
    id_max_abs_corr_fdr);

%%%%%%%%%%%%%%%%%%%%%%%%%%% BONFERRONI CORRECTED %%%%%%%%%%%%%%%%%%%%%%%%%%

% Inflation percentage of thresholded Bonferroni corrected correlations
type_m_errors.inflation_percentage_thresholded_bonferroni = ...
    aggregated_inflation_percentage_thresholded_bonferroni;

% Inflation percentage of Bonferroni corrected correlations
type_m_errors.inflation_percentage_bonferroni = ...
    aggregated_inflation_percentage_bonferroni;

% Average inflation percentage of thresholded Bonferroni corrected 
% correlations
type_m_errors.average_inflation_percentage_thresholded_bonferroni = ...
    average_inflation_thresholded_bonferroni;

% Average inflation percentage of Bonferroni corrected correlations
type_m_errors.average_inflation_percentage_bonferroni = ...
    average_inflation_bonferroni;

% Inflation percentage of the best correlation thresholded Bonferroni 
% corrected
type_m_errors.best_inflation_percentage_thresholded_bonferroni = ...
    aggregated_inflation_percentage_thresholded_bonferroni(:,...
    id_max_abs_corr_bonferroni);

% Inflation percentage of the best correlation Bonferroni corrected
type_m_errors.best_inflation_percentage_bonferroni = ...
    aggregated_inflation_percentage_bonferroni(:,...
    id_max_abs_corr_bonferroni);

% Timestamp
elapsed = seconds(toc(t0));
elapsed.Format = 'hh:mm:ss.SSS';
elapsed = string(elapsed);

% Provide feedback
fprintf('All done! Elapsed time: %s\n',...
    elapsed);

end
