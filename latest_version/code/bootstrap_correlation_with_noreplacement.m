function bootstrap_output = bootstrap_correlation_with_noreplacement(behavior, brain, alpha_level, sample_size_range, n_sample_sizes, n_bootstraps)

% behavior: a n-by-1 array storing simulated behavioral scores, where n
% represents the number of participants in the full sample. It can be
% generated using the generate_simulated_data.m function
%
% brain: a n-by-r matrix storing simulated brain data, where n
% represents the number of participants in the full sample and r the number
% of edges/ROIs/networks/components in the analysis. It can be generated
% using the generate_simulated_data.m function
%
% alpha_level: a number that defines the two-tails level of significance 
% (e.g., 0.05). This value is then adjusted for multiple comparisons based 
% on the number of edges/ROIs/components/networks using two popular p-value
% adjustment methods: Bonferroni and False Discovery Rate (Benjamini and
% Hochberg, 1995)
%
% sample_size_range: a 1-by-2 array defining the minimum and maximum sample
% sizes for which power is calculated. For instance [25 3000]. Of course
% the maximum number should be equal to or smaller than the number of
% observations in the full sample
%
% n_sample_sizes: the number of sample sizes in the sample_size_range for
% which power is computed. For instance 16
%
% n_bootstraps: number of resamplings at each sample size. This parameter 
% relates to the precision we have in computing power, type M, type S etc.
% The value reported in the article is 1000 and power is estimated with 
% 0.1% precision.

% coded by Luca Cecchetti & Giacomo Handjaras, IMT-Lucca, Italy
% vers 20220619

% Timestamp
t0 = tic;

% Begin of the analysis part
fprintf('\n----------------------------------\n')
fprintf('| Estimate Bootstrap Correlation |\n')
fprintf('----------------------------------\n')

% Number of tests
n_actual_correlations = size(brain,2);

% Number of participants in the full sample
n_full_sample = size(behavior,1);

% Number of steps in progress bar
step_progress = 10;

% Compute the correlation between brain and behavior in the full sample and
% estimate statistical significance
[r_full,p_full] = corr(brain,behavior);

% Get the correlation sign in the full sample
r_full_sign = sign(r_full);

% Correlations passing the uncorrected alpha threshold in the full sample
p_full_significant_uncorrected = p_full < alpha_level;

% Correlations passing the FDR corrected alpha threshold in the full sample
p_full_significant_fdr = fdr_bh(p_full,alpha_level,'pdep');

% Correlations passing the bonferroni corrected alpha threshold in the full
% sample
p_full_significant_bonferroni = ...
    p_full < (alpha_level/n_actual_correlations);

% Define an array storing the number of observations in each simulated
% sample size. We sample using the same parameters employed by Marek and
% Tervo-Clemmens. Thos used in the original article are [25, 33, 50, 70,
% 100, 135, 200, 265, 375, 525, 725, 1000, 1430, 2000, 2800]. If you would
% like to get something similar you can specify:
% sample_size_range = [25 3928] and n_sample_sizes = 16
sample_sizes = round(...
    logspace(...
    log10(sample_size_range(1)),...
    log10(sample_size_range(2)),...
    n_sample_sizes));

% Preallocate an array to store simulated brain-behavior correlations at
% different sample sizes and effect sizes and for each random resampling
simulated_corr = ...
    nan(n_bootstraps,n_sample_sizes,n_actual_correlations);

% Preallocate an array to store the level of significane for brain-behavior
% correlations at different sample sizes and effect sizes and for each
% random resampling
simulated_pvalue = ...
    nan(n_bootstraps,n_sample_sizes,n_actual_correlations);

% Preallocate an array to store brain-behavior correlations passing the FDR
% correction for multiple comparisons at different sample sizes and effect
% sizes and for each random resampling
simulated_significant_pvalue_fdr = ...
    false(n_bootstraps,n_sample_sizes,n_actual_correlations);

% For each random resampling
for b = 1:n_bootstraps
    
    % Randomize the order of simulated participants
    resampling_array = randperm(n_full_sample);
    
    % For each sample size
    for o = 1:n_sample_sizes
        
        % Number of observations in the selected sample size
        n_obs = sample_sizes(o);
        
        % Randomly select observations based on sample size. Please note we
        % are adding to the smallest group of participants more and more
        % individuals for each bootstrap based on the selected sample size.
        % This is because we don't want to introduce a random bias across
        % sample sizes pertaining to the same resampling
        random_sample = resampling_array(1:n_obs);
        
        % Select random simulated brain data
        brain_boot = brain(random_sample,:);
        
        % Select random simulated behavioral data
        behavior_boot = behavior(random_sample);
        
        % Estimate the correlation and its significance and store them into
        % separate matrices. Please note, here we computed statistical
        % significance based on parametric testing because of two reasons:
        % (1) data are generated from normal distributions, thus
        % assumptions for parametric testing are very likely to be met; (2)
        % if the assumptions are not met (because of any reason), we want
        % to be prone to false positives, as we are interested in getting
        % the number of small sample size studies (i.e., resamplings)
        % reporting an opposite effect as compared to the true one.
        [simulated_corr(b,o,:),simulated_pvalue(b,o,:)] = ...
            corr(brain_boot,behavior_boot);
        
        % Apply FDR correction for multiple comparisons
        simulated_significant_pvalue_fdr(b,o,:) = ...
            fdr_bh(squeeze(simulated_pvalue(b,o,:)),alpha_level,'pdep');
        
        
    end
    
    % Print feedback about progress
    if mod(b,floor(n_bootstraps/step_progress))==0
        
        % Time stamp
        elapsed = seconds(toc(t0));
        elapsed.Format = 'hh:mm:ss.SSS';
        elapsed = string(elapsed);
        
        % Provide feedback
        fprintf('Completed %d out of %d bootstraps - Elapsed time: %s\n',...
            b,n_bootstraps,elapsed);
        
    end
    
    
end

% Get the correlation sign at each sample size and random resampling
simulated_corr_sign = sign(simulated_corr);

% Threshold pvalues using the uncorrected alpha level
simulated_significant_pvalue_uncorrected = ...
    simulated_pvalue < alpha_level;    

% Apply Bonferroni correction for multiple comparisons
simulated_significant_pvalue_bonferroni = ...
    simulated_pvalue < (alpha_level/n_actual_correlations);

% Prepare output structure

% Sample sizes 
bootstrap_output.sample_sizes = sample_sizes;

% Correlations in the full sample
bootstrap_output.r_full = r_full;

% Sign of correlations in the full sample
bootstrap_output.r_full_sign = r_full_sign;

% Significance of the correlation in the full sample
bootstrap_output.p_full = p_full;

% Correlations passing uncorrected alpha level in the full sample
bootstrap_output.p_full_significant_uncorrected = ...
    p_full_significant_uncorrected;

% Correlations passing FDR corrected alpha level in the full sample
bootstrap_output.p_full_significant_fdr = ...
    p_full_significant_fdr;

% Correlations passing Bonferroni corrected alpha level in the full sample
bootstrap_output.p_full_significant_bonferroni = ...
    p_full_significant_bonferroni;

% Bootstrapped correlations at different sample sizes
bootstrap_output.simulated_corr = ...
    simulated_corr;

% Sign of bootstrapped correlations at different sample sizes
bootstrap_output.simulated_corr_sign = ...
    simulated_corr_sign;

% Significance of bootstrapped correlations at different sample sizes
bootstrap_output.simulated_pvalue = ...
    simulated_pvalue;

% Bootstrapped correlations passing uncorrected alpha threshold at 
% different sample sizes
bootstrap_output.simulated_significant_pvalue_uncorrected = ...
    simulated_significant_pvalue_uncorrected;

% Bootstrapped correlations passing FDR corrected alpha threshold at 
% different sample sizes
bootstrap_output.simulated_significant_pvalue_fdr = ...
    simulated_significant_pvalue_fdr;

% Bootstrapped correlations passing Bonferroni corrected alpha threshold at 
% different sample sizes
bootstrap_output.simulated_significant_pvalue_bonferroni = ...
    simulated_significant_pvalue_bonferroni;

% Timestamp
elapsed = seconds(toc(t0));
elapsed.Format = 'hh:mm:ss.SSS';
elapsed = string(elapsed);

% Provide feedback
fprintf('All done! Elapsed time: %s\n',...
    elapsed);



end
