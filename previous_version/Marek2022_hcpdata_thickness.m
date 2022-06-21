%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% coded by Luca Cecchetti & Giacomo Handjaras, IMT-Lucca, Italy
%%% vers 20220322
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc

addpath('bwas-main/');
addpath('additional_scripts/');

% Fix the seed for reproducibility
rng(14051983)

% Load HCP unrestricted data. We are not allowed to directly share HCP data
% but you can apply to get access to unrestricted data, which is sufficient
% to replicate the findings we report in the response letter. Change the
% name of the csv file accordingly
hcp_data = readtable('unrestricted_lucacecchetti_3_21_2022_1_1_51.csv');

% Get the index of columns storing Freesurfer cortical thickness values
idx_brain_data = startsWith(hcp_data.Properties.VariableNames,'FS_') & ...
    endsWith(hcp_data.Properties.VariableNames,'Thck');

% Get the ROI name of Freesurfer cortical thickness data
roiname_brain_data = hcp_data.Properties.VariableNames(...
    startsWith(hcp_data.Properties.VariableNames,'FS_') & ...
    endsWith(hcp_data.Properties.VariableNames,'Thck'));

% Set the significance level
alpha_level = 0.05;

% Number of random resamplings
n_boots = 1000;

% Number of sample sizes to be tested
n_sample_sizes = 16;

% Colormap for plotting results
plot_colormap = [189,201,225;...
    103,169,207;...
    2,129,138]./255;

% Get brain data from the HCP table
brain_data = table2array(hcp_data(:,idx_brain_data));

% Get behavioral data from HCP table. We will focus on
% CogCrystalComp_AgeAdj, similarly to Marek, Tervo-Clemmens et al.
behav_data = hcp_data.CogTotalComp_AgeAdj;

% Aggregate behavioral and brain data into a new table so that it will be
% easier to remove missing values
aggregated_table = cat(2,behav_data,brain_data);

% Remove missing values (participants)
aggregated_table_no_missing = rmmissing(aggregated_table);

% Get the number of participants
n_subs = size(aggregated_table_no_missing,1);

% Get the number of ROIs
n_rois = sum(idx_brain_data);

% Split behavioral data from brain data after missing values are removed
behav_data_no_missing = aggregated_table_no_missing(:,1);

% Split brain data from behavioral data after missing values are removed
brain_data_no_missing = aggregated_table_no_missing(:,2:end);

% Correlate brain and behavior in the full sample
[r_full, p_full] = corr(behav_data_no_missing,brain_data_no_missing);

% Find the strongest correlation (i.e., effect size) across all ROIs
[~,id_r_full_max] = max(abs(r_full));

% Provide feedback about the strongest relationship found
fprintf('The largest effect is %.5f\n',r_full(id_r_full_max));

% Find significant uncorrected correlations in the full sample
uncorrected_significant_p = p_full < alpha_level;

% Find significant fdr corrected correlations in the full sample (Benjamini
% and Hochberg method)
fdr_pdep_significant_p = fdr_bh(p_full, alpha_level, 'pdep');

% Find significant bonferroni corrected correlations in the full sample
bonferroni_significant_p = (p_full .* n_rois) < alpha_level;

% Array storing the sample sizes for the resampling procedure
sample_sizes = round(linspace(25,1000,n_sample_sizes));

% Preallocate an array to store correlations of resampled data
r_resampled = nan(n_boots,n_rois,n_sample_sizes);

% Preallocate an array to store raw p-values of resampled data
p_resampled = nan(n_boots,n_rois,n_sample_sizes);

% Preallocate an array to store significant uncorrected pvalues of
% resampled data
uncorrected_significant_p_resampled = false(n_boots,n_rois,n_sample_sizes);

% Preallocate an array to store significant fdr-corrected pvalues of
% resampled data
fdr_pdep_significant_p_resampled = false(n_boots,n_rois,n_sample_sizes);

% Preallocate an array to store significant bonferroni-corrected pvalues of
% resampled data
bonferroni_significant_p_resampled = false(n_boots,n_rois,n_sample_sizes);

% Variable to control verbosity of feedback: 10 means 10 steps
step_progress = 10;

% Time stamp
t0 = tic;

% For each random resampling
for i = 1:n_boots
    
    % Suffle the position of participants
    resampling_array = randperm(n_subs);
    
    % For each sample size
    for o = 1:n_sample_sizes
        
        % Get the number of observations in sample size
        n_obs_subsample = sample_sizes(o);
        
        % Get behavioral data after resampling
        behav_data_no_missing_resampled = ...
            behav_data_no_missing(resampling_array(1:n_obs_subsample));
        
        % Get brain data after resampling
        brain_data_no_missing_resampled = ...
            brain_data_no_missing(resampling_array(1:n_obs_subsample),:);
        
        % Correlate brain and behavior in resampled data
        [r_res, p_res] = corr(behav_data_no_missing_resampled,...
            brain_data_no_missing_resampled);
        
        % Save resampled correlation values
        r_resampled(i,:,o) = r_res;
        
        % Save raw pvalues for resampled correlation values
        p_resampled(i,:,o) = p_res;
        
        % Store logical array of pvalues passing (i.e., == 1) or not (i.e.,
        % == 0) the uncorrected statistical threhold
        uncorrected_significant_p_resampled(i,:,o) = ...
            p_res < alpha_level;
        
        % Store logical array of pvalues passing (i.e., == 1) or not (i.e.,
        % == 0) the fdr-corrected statistical threhold
        fdr_pdep_significant_p_resampled(i,:,o) = ...
            fdr_bh(p_res, alpha_level, 'pdep');
        
        % Store logical array of pvalues passing (i.e., == 1) or not (i.e.,
        % == 0) the bonferroni corrected statistical threhold
        bonferroni_significant_p_resampled(i,:,o) = ...
            (p_res .* n_rois) < alpha_level;
        
    end
    
    % Print feedback about progress
    if mod(i,floor(n_boots/step_progress))==0
        
        % Time stamp
        elapsed = seconds(toc(t0));
        elapsed.Format = 'hh:mm:ss.SSS';
        elapsed = string(elapsed);
        
        % Provide feedback
        fprintf('Completed %d out of %d bootstraps - Elapsed time: %s\n',...
            i,n_boots,elapsed);
        
    end
    
    
end

% Compute power for uncorrected pvalues. Power of a resampled version of
% the original data (i.e., at different n) is the number of times a
% significant ROI in the full sample is found to be significant in
% resampled data as well (i.e., smaller n group). Nota bene: this measure
% is ROI specific as each ROI correlates with behavior to some extent.
power_uncorrected = ...
    sum(uncorrected_significant_p_resampled(:,uncorrected_significant_p,:)) ...
    ./ n_boots .* 100;

% Compute power for fdr-corrected pvalues. Same as above but using fdr-bh
% correction method
power_fdr_pdep = ...
    sum(fdr_pdep_significant_p_resampled(:,fdr_pdep_significant_p,:)) ...
    ./ n_boots .* 100;

% Compute power for bonferroni-corrected pvalues. Same as above but using
% bonferroni correction method
power_bonferroni = ...
    sum(bonferroni_significant_p_resampled(:,bonferroni_significant_p,:)) ...
    ./ n_boots .* 100;

% Catenate power for each ROI and correction method
all_powers = cat(3,...
    squeeze(sum(uncorrected_significant_p_resampled)),...
    squeeze(sum(fdr_pdep_significant_p_resampled)),...
    squeeze(sum(bonferroni_significant_p_resampled)))./ n_boots .* 100;


% Plot power as a function of sample size and correction method for the
% strongest relationship between brain and behavior in the full sample.
figure;
plot(sample_sizes,squeeze(all_powers(id_r_full_max,:,1)),...
    'LineWidth',3.5,'Color',plot_colormap(1,:))
hold on
plot(sample_sizes,squeeze(all_powers(id_r_full_max,:,2)),...
    'LineWidth',3.5,'Color',plot_colormap(2,:))
plot(sample_sizes,squeeze(all_powers(id_r_full_max,:,3)),...
    'LineWidth',3.5,'Color',plot_colormap(3,:))
set(gca, 'XScale', 'linear',...
    'XTick',[300 600 900],...
    'XTickLabel',{'300','600','900'},...
    'YTick',0:20:100,...
    'TickDir','out',...
    'FontName','Arial',...
    'FontSize',22,...
    'YGrid','on');
xlim([min(sample_sizes),max(sample_sizes)])
ylim([0 100])
box off
set(gcf,'color',[1 1 1])
ylabel('Power (%)')
xlabel('Sample size')
legend({'p<0.05','q<0.05', 'p_B_o_n_f<0.05'},...
    'location','eastoutside','Box','off','FontSize',20)
axis square

% Uncomment the following line if you want to save to disk the figure in
% our reply to Marek, Tervo-Clemmens et al.
% For this to work you need a package that can be downloaded here:
% https://it.mathworks.com/matlabcentral/fileexchange/23629-export_fig
% 
% UNCOMMENT THIS:
% export_fig('hcp_best_roi_power_thickness_nonadjusted.png','-m6','-p0.05')

% Now we use the Marek, Tervo-Clemmens et al. function to compute type I,
% type II, type S and type M errors reported in Figure 3.
% The function is abcd_statisticalerrors.m and cen be downloaded from
% https://gitlab.com/DosenbachGreene/bwas
% Even though it is only seldom documented, it requires a variable storing
% correlation values for all resamplings, all sample sizes and all ROIs.
% Similarly it requires the pvalues associated to each resampled
% correlation. It also requires the correlation and pvalue computed in the
% full sample for each ROI. The last argument is the number of resamplings

% Correlation and pvalues of resampled data should be arranged so that the
% matrix is a regions-by-bootstraps-by-sample sizes
r_resampled_marek = permute(r_resampled,[2,1,3]);
p_resampled_marek = permute(p_resampled,[2,1,3]);

% Define the uncorrected and bonferroni corrected thresholds to be used in
% the power calculation according to Marek and colleagues
alpha_level_marek = [alpha_level/n_rois, alpha_level];

% Run power calculation as in Marek, Tervo-Clemmens et al. The power is
% computed for each ROI and then averaged across significant ROIs. Power is
% reported for uncorrected and bonferroni corrected pvalues.
[type1,type2,typem,types,types_pvals] = ...
    abcd_statisticalerrors_ch_mod(r_resampled_marek,...
    p_resampled_marek,...
    r_full',p_full',n_boots,alpha_level_marek);

% Power is defined as 100 minus type 2 error rates
power_marek = 100-type2;

% Plot power as a function of sample size and correction method for the
% strongest relationship between brain and behavior in the full sample.
figure;
plot(sample_sizes,power_marek(:,2),...
    'LineWidth',3.5,'Color',plot_colormap(1,:))
hold on
plot(sample_sizes,power_marek(:,1),...
    'LineWidth',3.5,'Color',plot_colormap(3,:))
set(gca, 'XScale', 'linear',...
    'XTick',[300 600 900],...
    'XTickLabel',{'300','600','900'},...
    'YTick',0:20:100,...
    'TickDir','out',...
    'FontName','Arial',...
    'FontSize',22,...
    'YGrid','on');
xlim([min(sample_sizes),max(sample_sizes)])
ylim([0 100])
box off
set(gcf,'color',[1 1 1])
ylabel('Power (%)')
xlabel('Sample size')
legend({'p<0.05','p_B_o_n_f<0.05'},...
    'location','eastoutside','Box','off','FontSize',20)
axis square

% Uncomment the following line if you want to save to disk the figure in
% our reply to Marek, Tervo-Clemmens et al.
% For this to work you need a package that can be downloaded here:
% https://it.mathworks.com/matlabcentral/fileexchange/23629-export_fig
% 
% UNCOMMENT THIS:
% export_fig('hcp_marek_power_thickness_nonadjusted.png','-m6','-p0.05')




