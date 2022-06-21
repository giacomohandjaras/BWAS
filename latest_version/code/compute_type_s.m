function type_s_errors = compute_type_s(bootstrap_output)
% coded by Luca Cecchetti & Giacomo Handjaras, IMT-Lucca, Italy
% vers 20220619

% Timestamp
t0 = tic;

% Begin of the analysis part
fprintf('\n------------------------------\n')
fprintf('| Estimate Type S Error Rate |\n')
fprintf('------------------------------\n')

% The best correlation in the full sample
[~, id_max_abs_corr] = max(abs(bootstrap_output.r_full));

% Number of bootstraps
n_bootstraps = size(bootstrap_output.simulated_corr_sign,1);

% Number of sample sizes
n_sample_sizes = size(bootstrap_output.simulated_corr_sign,2);

% Number of correlations
n_actual_correlations = size(bootstrap_output.simulated_corr_sign,3);

% Preallocate a matrix to store sign flips of unthresholded correlations
type_s_unthresholded_percentage = ...
    nan(n_sample_sizes,n_actual_correlations);

% Preallocate a matrix to store sign flips of thresholded uncorrected 
% correlations. With "thresholded correlations" we refer to those passing 
% significance level at each sample size and resampling
type_s_thresholded_uncorrected_percentage = ...
    nan(n_sample_sizes,n_actual_correlations);

% Preallocate a matrix to store sign flips of thresholded FDR corrected 
% correlations. With "thresholded correlations" we refer to those passing 
% significance level at each sample size and resampling
type_s_thresholded_fdr_percentage = ...
    nan(n_sample_sizes,n_actual_correlations);

% Preallocate a matrix to store sign flips of thresholded Bonferroni 
% corrected correlations. With "thresholded correlations" we refer to those 
% passing significance level at each sample size and resampling
type_s_thresholded_bonferroni_percentage = ...
    nan(n_sample_sizes,n_actual_correlations);

% For each sample size
for s = 1:n_sample_sizes
    
    % For each correlation
    for c = 1:n_actual_correlations
        
        % The sign of the correlation in the full sample is multiplied by
        % the sign of the correlation in each sample size and bootstrap. If
        % the output is 1 then signs are the same, if the sign is -1 then
        % signs are opposite
        temp_sign_flips = ...
            squeeze(bootstrap_output.simulated_corr_sign(:,s,c)) ...
            .* bootstrap_output.r_full_sign(c);
        
        % Get sign flip for correlations passing uncorrected significance
        % threshold at each sample size and resampling
        temp_sign_flips_uncorrected = ...
            temp_sign_flips(...
            bootstrap_output.simulated_significant_pvalue_uncorrected(...
            :,s,c));
        
        % Get sign flip for correlations passing FDR corrected significance
        % threshold at each sample size and resampling
        temp_sign_flips_fdr = ...
            temp_sign_flips(...
            bootstrap_output.simulated_significant_pvalue_fdr(:,s,c));
        
        % Get sign flip for correlations passing Bonferroni corrected 
        % significance threshold at each sample size and resampling
        temp_sign_flips_bonferroni = ...
            temp_sign_flips(...
            bootstrap_output.simulated_significant_pvalue_bonferroni(...
            :,s,c));                
        
        % Estimate the proportion of correlations with sign flip
        type_s_unthresholded_percentage(s,c) = ...
            sum(temp_sign_flips == -1) / n_bootstraps * 100;
        
        % Estimate the proportion of correlations passing uncorrected
        % threshold with sign flip
        type_s_thresholded_uncorrected_percentage(s,c) = ...
            sum(temp_sign_flips_uncorrected == -1) / n_bootstraps * 100;
        
        % Estimate the proportion of correlations passing FDR corrected
        % threshold with sign flip
        type_s_thresholded_fdr_percentage(s,c) = ...
            sum(temp_sign_flips_fdr == -1) / n_bootstraps * 100;
        
        % Estimate the proportion of correlations passing Bonferroni 
        % corrected threshold with sign flip
        type_s_thresholded_bonferroni_percentage(s,c) = ...
            sum(temp_sign_flips_bonferroni == -1) / n_bootstraps * 100;
                            
    end
            
end

% Sign flips for the best correlation at each sample size
best_sign_flip = ...
    squeeze(bootstrap_output.simulated_corr_sign(:,:,id_max_abs_corr)) ...
    .* bootstrap_output.r_full_sign(id_max_abs_corr);

% Proportion of sign flips for the best correlation
best_type_s_percentage = ...
    type_s_unthresholded_percentage(:,id_max_abs_corr);

% Proportion of sign flips for the best correlation considering resamplings
% passing the uncorrected threshold
best_type_s_thresholded_uncorrected_percentage = ...
    sum(best_sign_flip ...
    .* bootstrap_output.simulated_significant_pvalue_uncorrected...
    (:,:,id_max_abs_corr) == ...
    -1) ./ n_bootstraps .* 100;

% Proportion of sign flips for the best correlation considering resamplings
% passing the FDR corrected threshold
best_type_s_thresholded_fdr_percentage = ...
    sum(best_sign_flip ...
    .* bootstrap_output.simulated_significant_pvalue_fdr...
    (:,:,id_max_abs_corr) == ...
    -1) ./ n_bootstraps .* 100;

% Proportion of sign flips for the best correlation considering resamplings
% passing the Bonferroni corrected threshold
best_type_s_thresholded_bonferroni_percentage = ...
    sum(best_sign_flip ...
    .* bootstrap_output.simulated_significant_pvalue_bonferroni...
    (:,:,id_max_abs_corr) == ...
    -1) ./ n_bootstraps .* 100;

% Select only correlations passing uncorrected threshold in both full
% sample and each sample size
type_s_thresholded_uncorrected_percentage = ...
    type_s_thresholded_uncorrected_percentage...
    (:,bootstrap_output.p_full_significant_uncorrected);

% Select only correlations passing FDR corrected threshold in both full
% sample and each sample size
type_s_thresholded_fdr_percentage = ...
    type_s_thresholded_fdr_percentage...
    (:,bootstrap_output.p_full_significant_fdr);

% Select only correlations passing Bonferroni corrected threshold in both 
% full sample and each sample size
type_s_thresholded_bonferroni_percentage = ...
    type_s_thresholded_bonferroni_percentage...
    (:,bootstrap_output.p_full_significant_bonferroni);

% Mean type s error percentage for each sample size, considering
% correlations passing the uncorrected threshold in the full sample
average_type_s_unthresholded_percentage = ...
    mean(...
    type_s_unthresholded_percentage,2);

% Mean type s error percentage for each sample size, considering
% correlations passing the uncorrected threshold in the full sample
average_type_s_uncorrected_percentage = ...
    mean(...
    type_s_unthresholded_percentage...
    (:,bootstrap_output.p_full_significant_uncorrected),2);

% Mean type s error percentage for each sample size, considering
% correlations passing the fdr threshold in the full sample
average_type_s_fdr_percentage = ...
    mean(...
    type_s_unthresholded_percentage...
    (:,bootstrap_output.p_full_significant_fdr),2);

% Mean type s error percentage for each sample size, considering
% correlations passing the bonferroni threshold in the full sample
average_type_s_bonferroni_percentage = ...
    mean(...
    type_s_unthresholded_percentage...
    (:,bootstrap_output.p_full_significant_bonferroni),2);

% Mean type s error percentage for each sample size, considering
% correlations passing the uncorrected threshold in both the full sample
% and each resampling
average_type_s_thresholded_uncorrected_percentage = ...
    mean(type_s_thresholded_uncorrected_percentage,2);

% Mean type s error percentage for each sample size, considering
% correlations passing the FDR corrected threshold in both the full sample
% and each resampling
average_type_s_thresholded_fdr_percentage = ...
    mean(type_s_thresholded_fdr_percentage,2);

% Mean type s error percentage for each sample size, considering
% correlations passing the Bonferroni corrected threshold in both the full 
% sample and each resampling
average_type_s_thresholded_bonferroni_percentage = ...
    mean(type_s_thresholded_bonferroni_percentage,2);


% Prepare output

% Type S errors of the best correlation
type_s_errors.best_unthresholded = best_type_s_percentage;

% Type S errors of the best correlation for resamplings passing the
% uncorrected threshold
type_s_errors.best_thresholded_uncorrected = ...
    best_type_s_thresholded_uncorrected_percentage;

% Type S errors of the best correlation for resamplings passing the FDR
% corrected threshold
type_s_errors.best_thresholded_fdr = ...
    best_type_s_thresholded_fdr_percentage;

% Type S errors of the best correlation for resamplings passing the
% Bonferroni corrected threshold
type_s_errors.best_thresholded_bonferroni = ...
    best_type_s_thresholded_bonferroni_percentage;

% Type S errors of all correlations
type_s_errors.unthresholded = type_s_unthresholded_percentage;

% Average Type S error across all correlations
type_s_errors.average_unthresholded = ...
    average_type_s_unthresholded_percentage;

% Type S errors of correlations passing the uncorrected threshold
type_s_errors.uncorrected = type_s_unthresholded_percentage...
    (:,bootstrap_output.p_full_significant_uncorrected);

% Average Type S error across correlations passing the uncorrected 
% threshold
type_s_errors.average_uncorrected = average_type_s_uncorrected_percentage;

% Type S errors of correlations passing the FDR threshold
type_s_errors.fdr = type_s_unthresholded_percentage...
    (:,bootstrap_output.p_full_significant_fdr);

% Average Type S error across correlations passing the FDR threshold
type_s_errors.average_fdr = average_type_s_fdr_percentage;

% Type S errors of correlations passing the Bonferroni threshold
type_s_errors.bonferroni = type_s_unthresholded_percentage...
    (:,bootstrap_output.p_full_significant_bonferroni);

% Average Type S error across correlations passing the Bonferroni threshold
type_s_errors.average_bonferroni = average_type_s_bonferroni_percentage;

% Type S errors of correlations passing the uncorrected threshold in both
% full sample size and resamplings
type_s_errors.thresholded_uncorrected = ...
    type_s_thresholded_uncorrected_percentage;

% Average Type S errors across correlations passing the uncorrected 
% threshold in both full sample size and resamplings
type_s_errors.average_thresholded_uncorrected = ...
    average_type_s_thresholded_uncorrected_percentage;

% Type S errors of correlations passing the FDR corrected threshold in both
% full sample size and resamplings
type_s_errors.thresholded_fdr = ...
    type_s_thresholded_fdr_percentage;

% Average Type S errors across correlations passing the FDR corrected 
% threshold in both full sample size and resamplings
type_s_errors.average_thresholded_fdr = ...
    average_type_s_thresholded_fdr_percentage;

% Type S errors of correlations passing the Bonferroni corrected threshold 
% in both full sample size and resamplings
type_s_errors.thresholded_bonferroni = ...
    type_s_thresholded_bonferroni_percentage;

% Average Type S errors across correlations passing the Bonferroni 
% corrected threshold in both full sample size and resamplings
type_s_errors.average_thresholded_bonferroni = ...
    average_type_s_thresholded_bonferroni_percentage;

% Timestamp
elapsed = seconds(toc(t0));
elapsed.Format = 'hh:mm:ss.SSS';
elapsed = string(elapsed);

% Provide feedback
fprintf('All done! Elapsed time: %s\n',...
    elapsed);

end
