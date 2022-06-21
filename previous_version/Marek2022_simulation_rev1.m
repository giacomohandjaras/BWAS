%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% coded by Luca Cecchetti & Giacomo Handjaras, IMT-Lucca, Italy
%%% vers 20220323
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc

addpath('bwas-main/');
addpath('additional_scripts/');

% Fix the seed for reproducibility
rng(14051983)

% Determine whether results should be written to disk or not
save_results_to_disk = 'no';

% Number of ROIs in the Marek, Tervo-Clemmens et al. paper
n_rois_in_paper = 394;

% Define effect sizes (i.e., pearson's correlations) for the association
% between brain and behavior. The idea is to create an array of values
% comparable to those of figure 1A: relationship between cognition and
% cortical thickness. Thus, we selected effect sizes from a normal
% distribution with 0 mean but 0.0336 standard deviation, which we
% estimated using source data provided by Nature. To prove the difference
% between ROI-based power calculation and the average method proposed by
% Marek, Tervo-Clemmens et al., we simulate a very strong correlation
% between brain and behavior (i.e., r = 0.53) detectable with 80% power
% with n = 25 and compare the two procedures.

% Generate the random correlation values with a distribution comparable to
% figure 1A. The idea is that each value is a ROI
ground_truth_correlations = 0 + 0.0336 .* randn(1,n_rois_in_paper);

% Position of the ROI of interest in the ROIs array
ROI_of_interest = 1;

% Replace the first element of the random correlation with a very strong
% effect size, detectable using n = 25 with 80% of power. In other terms,
% the behavior has 0.53 correlation with ROI 1 and much smaller correlation
% with all other ROIs
ground_truth_correlations(ROI_of_interest) = 0.5300;

% Compute the number tested ROIs
n_ground_truth_correlations = numel(ground_truth_correlations);

% Define the number of comparisons based on the Marek, Tervo-Clemmens et
% al. ROI approach. The value 1 is used for computing significance using
% uncorrected pvalues
comparisons = [1, n_rois_in_paper];

% Number of analysis levels
n_comparisons = numel(comparisons);

% The size of the population
n_population = 100000;

% Overall number of participants (i.e., full sample size)
n_subs = 3928;

% Number of steps in progress bar
step_progress = 10;

% Noise level
sigma = .1;

% Preallocate covariance matrix of independent variables
covariance_structure = eye(n_rois_in_paper + 1);

% Impose covariance between first variable and all other variables based
% on ground truth correlations.
covariance_structure(2:end,1) = ground_truth_correlations;
covariance_structure(1,2:end) = ground_truth_correlations;

% Add noise to the covariance matrix
sigma_dep = sigma.^2 .* covariance_structure;

% Generate random variables using a multivariate normal distribution
data = mvnrnd(zeros(n_rois_in_paper+1,1),sigma_dep,n_population);

% Select the full sample size
data = data(1:n_subs,:);

% First column is random behavioral scores: each element simulates the
% behavioral performance of an individual in a test
behavior = data(:,1);

% All other columns are brain data (i.e., cortical thickness of multiple
% participants - rows, by regions of interest - columns)
brain = data(:,2:end);

% Compute the correlation and statistical significance between brain and
% behavior for the full sample
[r_full,p_full] = corr(brain,behavior);

% Define an array storing the number of observations in each simulated
% sample size. We sample using the same parameters employed by Marek and
% Tervo-Clemmens.
sample_sizes = [25, 33, 50, 70, 100, 135, 200, 265,...
    375, 525, 725, 1000, 1430, 2000, 2800, 3604];

% Number of sample sizes
n_sample_sizes = numel(sample_sizes);

% Number of resampling. The same value used in the paper
n_bootstraps = 1000;

% Define an array to store alpha levels for significance testing. These are
% the classic *uncorrected* thresholds. Bonferroni correction will be
% applied to each of the thresholds based on the number of comparisons in
% the corresponding analysis level array
alpha_levels = [0.05 0.01];

% Number of alpha levels
n_alpha_levels = numel(alpha_levels);

% Define confidence interval range 1 [min max]
prctile_range_1 = [2.5 97.5];

% Define confidence interval range 2 [min max]
prctile_range_2 = [0.5 99.5];

% Preallocate an array to store simulated brain-behavior correlations at
% different sample sizes and effect sizes and for each random resampling
simulated_corr = ...
    nan(n_bootstraps,n_sample_sizes,n_ground_truth_correlations);

% Preallocate an array to store the level of significane for brain-behavior
% correlations at different sample sizes and effect sizes and for each
% random resampling
simulated_pvalue = ...
    nan(n_bootstraps,n_sample_sizes,n_ground_truth_correlations);

% Begin of the analysis part
fprintf('----------------------\n')
fprintf('| Running Simulation |\n')
fprintf('----------------------\n')

% Time stamp
t0 = tic;

% For each random resampling
for b = 1:n_bootstraps
    
    % Randomize the order of simulated participants
    resampling_array = randperm(n_subs);
    
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

% Get the sign of the correlation coefficient for each bootstrap, sample
% size and tested effect
sign_simulated_corr = sign(simulated_corr);

% Get the proportion of correlations with opposite sign with respect to
% ground truth as a function of effect size and sample size
percentage_correlation_opposite_sign = ...
    (sum(sign_simulated_corr<0,1))./n_bootstraps.*100;

% Get the proportion of correlations with same sign with respect to
% ground truth as a function of effect size and sample size
percentage_correlation_same_sign = ...
    (sum(sign_simulated_corr>0,1))./n_bootstraps.*100;

% Preallocate a matrix to store significant pvalues at different alpha
% levels for correlations of opposite sign as compared to the ground truth
percentage_significant_correlation_opposite_sign = ...
    nan(n_sample_sizes,n_ground_truth_correlations,...
    n_alpha_levels,n_comparisons);

% Preallocate a matrix to store significant pvalues at different alpha
% levels for correlations with the same sign of the ground truth
percentage_significant_correlation_same_sign = ...
    nan(n_sample_sizes,n_ground_truth_correlations,...
    n_alpha_levels,n_comparisons);

% For each alpha level
for p = 1:n_alpha_levels
    
    % For each level of analysis (i.e., uncorrected, network-based,
    % ROI-based and vertex-based)
    for c = 1:n_comparisons
        
        % Estimate the (un)corrected alpha threshold
        p_thr = alpha_levels(p)/comparisons(c);
        
        % Get the matrix of significant pvalues, expressing also the
        % correlation sign: elements == 0 are not significant, elements ==
        % 1 are significant and have the same sign of the ground truth,
        % elements == -1 are significant and opposite sign as compared to
        % the ground truth
        significant_p_w_corr_sign = ...
            sign_simulated_corr .* (simulated_pvalue < p_thr);
        
        % Get the proportion of significant correlations with opposite sign
        % with respect to ground truth as a function of effect size and
        % sample size
        percentage_significant_correlation_opposite_sign(:,:,p,c) = ...
            squeeze(sum(significant_p_w_corr_sign<0,1))./n_bootstraps.*100;
        
        % Get the proportion of significant correlations with same sign of
        % the ground truth as a function of effect size and sample size
        percentage_significant_correlation_same_sign(:,:,p,c) = ...
            squeeze(sum(significant_p_w_corr_sign>0,1))./n_bootstraps.*100;
        
    end
    
end

% Store results to disk
if startsWith(save_results_to_disk,'y')
    
    save('marek_tervoclemmens_simulation_results.mat','-v7.3')
    
end

% Plot the power as a function of sample size and correction method for the
% simulated strongest association between brain and behavior (i.e., ROI 1).
% We expect 80% power using the uncorrected p<0.05 and n=25.
plot_alpha_level = 1;
plot_ROI = ROI_of_interest;

% Colormap for plotting results
plot_colormap = [189,201,225;...
    103,169,207;...
    2,129,138]./255;


% Plot power as a function of sample size in correlations having the same
% sign of the ground truth across resamplings. For large effect sizes this
% is the same as computing the power over all significant correlations, as
% the possibility to get sign flip is non-existent
figure;
P1 = plot(sample_sizes,...
    percentage_significant_correlation_same_sign(:,...
    plot_ROI,plot_alpha_level,1),...
    'LineWidth',3.5,'Color',plot_colormap(1,:));
hold on
P2 = plot(sample_sizes,...
    percentage_significant_correlation_same_sign(:,...
    plot_ROI,plot_alpha_level,2),...
    'LineWidth',3.5,'Color',plot_colormap(3,:));
set(gca, 'XScale', 'log',...
    'XTick',[25 200 2000],...
    'XTickLabel',{'25','200','2000'},...
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
legend([P1,P2],...
    {'p<0.05','p_B_o_n_f<0.05'},'location','eastoutside','Box','off',...
    'FontSize',20)

% Uncomment the following line if you want to save to disk the figure in
% our reply to Marek, Tervo-Clemmens et al.
% For this to work you need a package that can be downloaded here:
% https://it.mathworks.com/matlabcentral/fileexchange/23629-export_fig
%
% UNCOMMENT THIS:
% filename_figure = ...
%     sprintf('simulation_best_roi_power_r_%.2f_n_%d.png',...
%     ground_truth_correlations(ROI_of_interest),n_subs);
% 
% export_fig(filename_figure,'-m6','-p0.05')


% Plot power as a function of sample size in correlations having opposite
% sign with respect to the ground truth across resamplings. This is also
% known as type III or type S error. For large effect sizes this is always
% zero. The idea is that if you set ground_truth_correlations to, for
% instance, 0.16 you will get a simulation for the largest effect size
% reported in the BWAS paper.
figure;
P3 = plot(sample_sizes,...
    percentage_significant_correlation_opposite_sign(:,...
    plot_ROI,plot_alpha_level,1),...
    'LineWidth',3.5,'Color',plot_colormap(1,:));
hold on
P4 = plot(sample_sizes,...
    percentage_significant_correlation_opposite_sign(:,...
    plot_ROI,plot_alpha_level,2),...
    'LineWidth',3.5,'Color',plot_colormap(3,:));
set(gca, 'XScale', 'log',...
    'XTick',[25 200 2000],...
    'XTickLabel',{'25','200','2000'},...
    'YTick',0:20:100,...
    'TickDir','out',...
    'FontName','Arial',...
    'FontSize',22,...
    'YGrid','on');
xlim([min(sample_sizes),max(sample_sizes)])
ylim([0 100])
box off
set(gcf,'color',[1 1 1])
ylabel('Type S Error (%)')
xlabel('Sample size')
legend([P3,P4],...
    {'p<0.05','p_B_o_n_f<0.05'},'location','eastoutside','Box','off',...
    'FontSize',20)

% Uncomment the following line if you want to save to disk the figure in
% our reply to Marek, Tervo-Clemmens et al.
% For this to work you need a package that can be downloaded here:
% https://it.mathworks.com/matlabcentral/fileexchange/23629-export_fig
%
% UNCOMMENT THIS:
% filename_figure = ...
%     sprintf('simulation_best_roi_typeS_r_%.2f_n_%d.png',...
%     ground_truth_correlations(ROI_of_interest),n_subs);
% 
% export_fig(filename_figure,'-m6','-p0.05')


% Colormap for plotting sign results
plot_colormap_sign = [251,128,114;...
    128,177,211]./255;

% Plot of the number of correlations with same vs opposite sign as compared
% to the correlation estimated in the full sample, across resamplings. For
% large effect sizes the lines should be maximally separated, whereas for
% negligible effects they should meet at the center of the y axis (i.e.,
% 50%)
figure;
C = plot(sample_sizes,...
    percentage_correlation_same_sign(:,:,plot_ROI),...
    'LineWidth',3.5,'Color',plot_colormap_sign(1,:));
hold on
O = plot(sample_sizes,...
    percentage_correlation_opposite_sign(:,:,plot_ROI),...
    'LineWidth',3.5,'Color',plot_colormap_sign(2,:));
set(gca, 'XScale', 'log',...
    'XTick',[25 200 2000],...
    'XTickLabel',{'25','200','2000'},...
    'YTick',0:20:100,...
    'TickDir','out',...
    'FontName','Arial',...
    'FontSize',22,...
    'YGrid','on');
xlim([min(sample_sizes),max(sample_sizes)])
ylim([0 100])
box off
set(gcf,'color',[1 1 1])
ylabel('Percentage of Correlations')
xlabel('Sample size')
legend([C,O],...
    {'Correct Sign','Opposite Sign'},'location','eastoutside','Box','off',...
    'FontSize',20)
% Uncomment the following line if you want to save to disk the figure in
% our reply to Marek, Tervo-Clemmens et al.
% For this to work you need a package that can be downloaded here:
% https://it.mathworks.com/matlabcentral/fileexchange/23629-export_fig
%
% UNCOMMENT THIS:
% filename_figure = ...
%     sprintf('simulation_best_roi_signflip_r_%.2f_n_%d.png',...
%     ground_truth_correlations(ROI_of_interest),n_subs);
% 
% export_fig(filename_figure,'-m6','-p0.05')

% Select correlation values for the ROI of interest
selected_simulated_corr = simulated_corr(:,:,plot_ROI);

% Select correlation values for the ROI of interest
selected_simulated_pvalue = simulated_pvalue(:,:,plot_ROI);

% Compute the average effect across random resamplings and sample sizes
avg_observed_correlations = squeeze(mean(selected_simulated_corr,1));

% Compute the lower boundary of CI 1 of the ROI across random
% resamplings and sample sizes
lower_prctile_1_observed_correlations = ...
    squeeze(prctile(selected_simulated_corr,prctile_range_1(1),1));

% Compute the upper boundary of CI 1 of the ROI across random
% resamplings and sample sizes
upper_prctile_1_observed_correlations = ...
    squeeze(prctile(selected_simulated_corr,prctile_range_1(2),1));

% Compute the lower boundary of CI 2 of the ROI across random
% resamplings and sample sizes
lower_prctile_2_observed_correlations = ...
    squeeze(prctile(selected_simulated_corr,prctile_range_2(1),1));

% Compute the upper boundary of CI 2 of the ROI across random
% resamplings and sample sizes
upper_prctile_2_observed_correlations = ...
    squeeze(prctile(selected_simulated_corr,prctile_range_2(2),1));

% Compute the max of the ROI across random resamplings and sample sizes
max_observed_correlations = ...
    squeeze(max(selected_simulated_corr,[],1));

% Compute the min of the ROI across random resamplings and sample sizes
min_observed_correlations = ...
    squeeze(min(selected_simulated_corr,[],1));

% x values for plotting max-min as a shaded area
shaded_x=[sample_sizes,fliplr(sample_sizes)];

% y values for plotting max-min as a shaded area
shaded_y=[max_observed_correlations,fliplr(min_observed_correlations)];

% Colormap for replicating Figure 1 with simulated data
plot_colormap_simulate = [106,106,106]./255;

% Replicating Marek et al., 2022 Figure 1 with simulated data
figure;
S = fill(shaded_x,shaded_y,...
    plot_colormap_simulate,...
    'FaceAlpha',0.3,...
    'EdgeColor','none');
hold on
A = plot(sample_sizes, avg_observed_correlations,...
    'LineStyle','-',...
    'LineWidth',1.2,...
    'Color',plot_colormap_simulate);
L95 = plot(sample_sizes, lower_prctile_1_observed_correlations,...
    'LineStyle',':',...
    'LineWidth',1,...
    'Color',[.2 .2 .2]);
U95 = plot(sample_sizes, upper_prctile_1_observed_correlations,...
    'LineStyle',':',...
    'LineWidth',1,...
    'Color',[.2 .2 .2]);
L99 = plot(sample_sizes, lower_prctile_2_observed_correlations,...
    'LineStyle','--',...
    'LineWidth',1,...
    'Color',[.2 .2 .2]);
U99 = plot(sample_sizes, upper_prctile_2_observed_correlations,...
    'LineStyle','--',...
    'LineWidth',1,...
    'Color',[.2 .2 .2]);
set(gca, 'XScale', 'log',...
    'XTick',[25 200 2000],...
    'XTickLabel',{'25','200','2000'},...
    'YTick',-1:0.5:1,...
    'TickDir','out',...
    'FontName','Arial',...
    'FontSize',22,...
    'YGrid','on');
xlim([min(sample_sizes),max(sample_sizes)])
ylim([-1 1])
box off
set(gcf,'color',[1 1 1])
ylabel('Correlation (r)')
xlabel('Sample size')
legend([S,A,L95,L99],...
    {'Min-Max','Avg','95CI','99CI'},'location','eastoutside','Box','off',...
    'FontSize',20)
% Uncomment the following line if you want to save to disk the figure in
% our reply to Marek, Tervo-Clemmens et al.
% For this to work you need a package that can be downloaded here:
% https://it.mathworks.com/matlabcentral/fileexchange/23629-export_fig
%
% UNCOMMENT THIS:
% filename_figure = ...
%     sprintf('simulation_best_roi_effectsize_r_%.2f_n_%d.png',...
%     ground_truth_correlations(ROI_of_interest),n_subs);
% 
% export_fig(filename_figure,'-m6','-p0.05')


% Preallocate a matrix to store the proportion of significant positive and
% negative associations according to the alpha level and correction for
% multiple comparisons. Please note the true effect (i.e., true_r) is
% positive, thus we should not observe any significant pvalue for negative
% (spurious) correlations. This is done for the ROI of interest
% First dimension: proportion of positive (i.e., 1st row) vs negative
% (i.e., 2nd row) significant correlations
% Second dimension: alpha level
% Third dimension: correction method
% Fourth dimension: sample sizes
sign_p = nan(2,...
    n_alpha_levels,n_comparisons,...
    numel(sample_sizes));

% Store power and type S error at different sample sizes and levels of
% significance in a cell array. Each element in the cell is a table.
% This should be human readable
sign_p_cell = cell(n_alpha_levels*n_comparisons,1);

% An array to store the number of negative (i.e., sign-flipped)
% correlations.
n_neg_corr_all = nan(1,numel(sample_sizes));

% For each sample size
for o = 1:numel(sample_sizes)
    
    % Identify positive correlations
    pos_corr = selected_simulated_corr(:,o)>0;
    
    % Identify negative correlations
    neg_corr = selected_simulated_corr(:,o)<0;
    
    % Count how many positive correlations
    n_pos_corr = sum(pos_corr);
    
    % Count how many negative correlations
    n_neg_corr = sum(neg_corr);
    
    % Store the number of negative correlations
    n_neg_corr_all(o) = n_neg_corr;
    
    % For each alpha level
    for p = 1:n_alpha_levels
        
        % For each level of analysis (i.e., uncorrected, ROI-based)
        for c = 1:n_comparisons
            
            % Select pvalues passing the significance threshold for
            % positive correlations
            pos_corr_sign = selected_simulated_pvalue(pos_corr,o) ...
                < (alpha_levels(p) / comparisons(c));
            
            % Select pvalues passing the significance threshold for
            % negative correlations
            neg_corr_sign = selected_simulated_pvalue(neg_corr,o) ...
                < (alpha_levels(p) / comparisons(c));
            
            % Number of significant pvalues for positive correlations
            n_pos_corr_sign_p = sum(pos_corr_sign);
            
            % Number of significant pvalues for negative correlations
            n_neg_corr_sign_p = sum(neg_corr_sign);
            
            % Store positive significant correlations
            sign_p(1,p,c,o) = n_pos_corr_sign_p;
            
            % Store negative significant correlations
            sign_p(2,p,c,o) = n_neg_corr_sign_p;
            
        end
        
    end
    
end

% Compute the power and type S error rate as percentages
sign_p = sign_p./n_bootstraps.*100;

% Variable to store tables in different cell array (this part is ugly as
% shit)
t=1;

% For each alpha level
for p = 1:n_alpha_levels
    
    % For each multiple comparisons type
    for c = 1:n_comparisons
        
        % Get the adjusted threshold
        p_thr = alpha_levels(p)/comparisons(c);
        
        % Convert the matrix in a table and add headers
        sign_p_cell{t} = ...
            array2table(squeeze(sign_p(:,p,c,:)),...
            'VariableNames',...
            strcat('N = ',{' '},strsplit(num2str(sample_sizes))),...
            'RowNames',{sprintf('Power (%%) - alpha level: %.5f',p_thr),...
            sprintf('Type S (%%) - alpha level: %.5f',p_thr)});
        
        % Increase the counter
        t = t+1;
        
    end
    
end

% Display tables storing power and type S error at all tested thresholds
% and sample sizes for the region of interest
sign_p_cell{:}

%%

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
r_resampled_marek = permute(simulated_corr,[3,1,2]);
p_resampled_marek = permute(simulated_pvalue,[3,1,2]);

% Define the uncorrected and bonferroni corrected thresholds to be used in
% the power calculation according to Marek and colleagues
alpha_level_marek = [alpha_levels(plot_alpha_level)/n_rois_in_paper,...
    alpha_levels(plot_alpha_level)];

% Run power calculation as in Marek, Tervo-Clemmens et al. The power is
% computed for each ROI and then averaged across significant ROIs. Power is
% reported for uncorrected and bonferroni corrected pvalues.
[type1,type2,typem,types,types_pvals] = ...
    abcd_statisticalerrors_ch_mod(r_resampled_marek,...
    p_resampled_marek,...
    r_full,p_full,n_bootstraps,alpha_level_marek);

% Power is defined as 100 minus type 2 error rates
power_marek = 100-type2;

% Plot power as a function of sample size and correction method for the
% strongest relationship between brain and behavior in the full sample.
figure;
P5 = plot(sample_sizes,power_marek(:,2),...
    'LineWidth',3.5,'Color',plot_colormap(1,:));
hold on
P6 = plot(sample_sizes,power_marek(:,1),...
    'LineWidth',3.5,'Color',plot_colormap(3,:));
set(gca, 'XScale', 'log',...
    'XTick',[25 200 2000],...
    'XTickLabel',{'25','200','2000'},...
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
legend([P5,P6],...
    {'p<0.05','p_B_o_n_f<0.05'},'location','eastoutside','Box','off',...
    'FontSize',20)

% Uncomment the following line if you want to save to disk the figure in
% our reply to Marek, Tervo-Clemmens et al.
% For this to work you need a package that can be downloaded here:
% https://it.mathworks.com/matlabcentral/fileexchange/23629-export_fig
%
% UNCOMMENT THIS:
% filename_figure = ...
%     sprintf('simulation_marek_power_r_%.2f_n_%d.png',...
%     ground_truth_correlations(ROI_of_interest),n_subs);
% 
% export_fig(filename_figure,'-m6','-p0.05')



