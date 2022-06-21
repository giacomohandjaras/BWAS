function statistical_power = compute_power(bootstrap_output)
% coded by Luca Cecchetti & Giacomo Handjaras, IMT-Lucca, Italy
% vers 20220619

% Timestamp
t0 = tic;

% Begin of the analysis part
fprintf('\n------------------\n')
fprintf('| Estimate Power |\n')
fprintf('------------------\n')

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

% Restrict the computation of power to the best correlation in the full
% sample passing the uncorrected threshold
best_power_uncorrected = ...
    power_uncorrected(:,id_max_abs_corr);

% Restrict the computation of power to correlations passing the significant
% uncorrected threshold in the full sample
power_uncorrected = ...
    power_uncorrected(:,bootstrap_output.p_full_significant_uncorrected);

% Compute average power for correlations passing the uncorrected
% significance threshold
average_power_uncorrected = mean(power_uncorrected,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FDR CORRECTED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute power at multiple sample sizes for FDR corrected significant 
% correlations
power_fdr = ...
    squeeze(...
    sum(bootstrap_output.simulated_significant_pvalue_fdr) ./ ...
    n_bootstraps .* 100);

% Restrict the computation of power to the best correlation in the full
% sample passing the FDR corrected threshold
best_power_fdr = ...
    power_fdr(:,id_max_abs_corr);

% Restrict the computation of power to correlations passing the FDR 
% significance threshold in the full sample
power_fdr = ...
    power_fdr(:,bootstrap_output.p_full_significant_fdr);

% Compute average power for correlations passing the FDR corrected
% significance threshold
average_power_fdr = mean(power_fdr,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%% BONFERRONI CORRECTED %%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute power at multiple sample sizes for Bonferroni corrected 
% significant correlations
power_bonferroni = ...
    squeeze(...
    sum(bootstrap_output.simulated_significant_pvalue_bonferroni) ./ ...
    n_bootstraps .* 100);

% Restrict the computation of power to the best correlation in the full
% sample passing the Bonferroni corrected threshold
best_power_bonferroni = ...
    power_bonferroni(:,id_max_abs_corr);

% Restrict the computation of power to correlations passing the bonferroni 
% significance threshold in the full sample
power_bonferroni = ...
    power_bonferroni(:,bootstrap_output.p_full_significant_bonferroni);

% Compute average power for correlations passing the Bonferroni corrected
% significance threshold
average_power_bonferroni = mean(power_bonferroni,2);

% Prepare output

% Power at multiple sample sizes for uncorrected significant correlations
% in the full sample
statistical_power.uncorrected = power_uncorrected;

% Power for the best correlation passing the uncorrected threshold
statistical_power.best_uncorrected = best_power_uncorrected;

% Average power for uncorrected correlations
statistical_power.average_uncorrected = average_power_uncorrected;

% Power at multiple sample sizes for uncorrected significant correlations
% in the full sample
statistical_power.fdr = power_fdr;

% Power for the best correlation passing the FDR corrected threshold
statistical_power.best_fdr = best_power_fdr;

% Average power for FDR corrected correlations
statistical_power.average_fdr = average_power_fdr;

% Power at multiple sample sizes for uncorrected significant correlations
% in the full sample
statistical_power.bonferroni = power_bonferroni;

% Power for the best correlation passing the Bonferroni corrected threshold
statistical_power.best_bonferroni = best_power_bonferroni;

% Average power for Bonferroni corrected correlations
statistical_power.average_bonferroni = average_power_bonferroni;

% Timestamp
elapsed = seconds(toc(t0));
elapsed.Format = 'hh:mm:ss.SSS';
elapsed = string(elapsed);

% Provide feedback
fprintf('All done! Elapsed time: %s\n',...
    elapsed);

end
