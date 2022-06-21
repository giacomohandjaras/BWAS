function false_positive_rate = compute_fpr(bootstrap_output)
% coded by Luca Cecchetti & Giacomo Handjaras, IMT-Lucca, Italy
% vers 20220619

% Timestamp
t0 = tic;

% Begin of the analysis part
fprintf('\n----------------------------\n')
fprintf('| Estimate False Positives |\n')
fprintf('----------------------------\n')

% Number of bootstraps
n_bootstraps = size(bootstrap_output.simulated_corr_sign,1);

% The best correlation in the full sample
[~, id_max_abs_corr] = max(abs(bootstrap_output.r_full));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UNCORRECTED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute power at multiple sample sizes for uncorrected significant 
% correlations
power_uncorrected = ...
    squeeze(...
    sum(bootstrap_output.simulated_significant_pvalue_uncorrected) ./ ...
    n_bootstraps .* 100);

% Restrict the computation of power to correlations not passing the 
% threshold for uncorrected significance in the full sample
fpr_uncorrected = ...
    power_uncorrected(:,~bootstrap_output.p_full_significant_uncorrected);

% Compute average false positives for correlations passing the uncorrected
% significance threshold
average_fpr_uncorrected = mean(fpr_uncorrected,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FDR CORRECTED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute power at multiple sample sizes for FDR corrected significant 
% correlations
power_fdr = ...
    squeeze(...
    sum(bootstrap_output.simulated_significant_pvalue_fdr) ./ ...
    n_bootstraps .* 100);

% Restrict the computation of power to correlations not passing the 
% threshold for FDR corrected significance in the full sample
fpr_fdr = ...
    power_fdr(:,~bootstrap_output.p_full_significant_fdr);

% Compute average false positives for correlations passing the FDR
% corrected significance threshold
average_fpr_fdr = mean(fpr_fdr,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%% BONFERRONI CORRECTED %%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute power at multiple sample sizes for Bonferroni corrected 
% significant correlations
power_bonferroni = ...
    squeeze(...
    sum(bootstrap_output.simulated_significant_pvalue_bonferroni) ./ ...
    n_bootstraps .* 100);

% Restrict the computation of power to correlations not passing the 
% threshold for Bonferroni corrected significance in the full sample
fpr_bonferroni = ...
    power_bonferroni(:,~bootstrap_output.p_full_significant_bonferroni);

% Compute average false positives for correlations passing the Bonferroni
% corrected significance threshold
average_fpr_bonferroni = mean(fpr_bonferroni,2);

% Prepare output

% False positives rate at multiple sample sizes for uncorrected significant
% correlations in the full sample
false_positive_rate.uncorrected = fpr_uncorrected;

% Average power for uncorrected correlations
false_positive_rate.average_uncorrected = average_fpr_uncorrected;

% False positives rate at multiple sample sizes for FDR corrected 
% significant correlations in the full sample
false_positive_rate.fdr = fpr_fdr;

% Average power for FDR corrected correlations
false_positive_rate.average_fdr = average_fpr_fdr;

% False positives rate at multiple sample sizes for Bonferroni corrected 
% significant correlations in the full sample
false_positive_rate.bonferroni = fpr_bonferroni;

% Average power for Bonferroni corrected correlations
false_positive_rate.average_bonferroni = average_fpr_bonferroni;

% Timestamp
elapsed = seconds(toc(t0));
elapsed.Format = 'hh:mm:ss.SSS';
elapsed = string(elapsed);

% Provide feedback
fprintf('All done! Elapsed time: %s\n',...
    elapsed);

end
