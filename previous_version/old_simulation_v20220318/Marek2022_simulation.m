%% Are small sample-size MRI studies useful? A simulated quantitative evaluation
%
% Authors:
% Luca Cecchetti, Social and Affective Neuroscience Group, IMT School for
% Advanced Studies Lucca, Italy
% Giacomo Handjaras, Social and Affective Neuroscience Group, IMT School
% for Advanced Studies Lucca, Italy 
%
% In this script we simulate data based on Figure 1 of Marek et al. (2022)
% "Reproducible brain-wide association studies require thousands of
% individuals". In Figure 1F, the authors show that if a true positive
% association between a behavioral measure and brain function/structure
% exists (largest r ~ +0.16 for the cognitive ability to DMN connectivity
% correlation) it is likely that small sample size studies will detect an
% effect with the opposite sign (i.e., negative correlation for the
% abovementioned relationship) because of sampling error. This, obviously, 
% would have dramatic consequences in terms of interpretation of findings,
% questioning the utility of small sample size brain-behavior association 
% studies.
% Nonetheless, what is not clear (IOHO not even from figure 3c) is the
% probability of observing a statistically significant opposite sign
% relationship in small sample-size studies. This point is crucial in the
% debate triggered by Marek et al., paper on the utility of small n
% studies, because the predicted effect of being underpowered is to fail in 
% detecting a true effect, rather than reporting a false positive effect 
% having opposite direction with respect to the true relationship.
% We try to address this issue with the following simulation.
% Note: we also used real data used to generate figure 1F, provided by
% Nature as "Source Data" (see Marke2022_actualdata.m script).
%
% RESULTS AND TAKE HOME: In simulated data (comparable to Figure 1 of Marek
% et al.), we do observe that sign flips in correlation coefficients are 
% likely to occurr when the sample size is small (e.g., 25, 33, 50). 
% However, across random resamplings, very few of these reverse 
% correlations pass the statistical threshold for uncorrected significance 
% (1 or 2 out of 1000 when true effect is r = 0.2) and none the bonferroni 
% corrected threshold (for the same effect size). Even when considering 
% top 1% effect sizes (e.g., >=|0.06|) the probability of observing 
% significant reverse correlations is extremely small (i.e., 1.4% of times 
% out of 1000 resamplings) for the uncorrected threshold and null for the 
% bonferroni corrected one. Therefore, small sample size brain-behavior 
% association studies are more likely to fail in detecting a true effect,
% or overestimating the effect size of a true relationship, rather than
% detecting spurious opposite sign correlations. This problem, of course 
% relates to statistical power of small n researches, which hampers any
% field of science.

clear
close all
clc

% Fix the seed for reproducibility
rng(14051983)

% Actual RSFC - cognitive ability correlation in Figure 1F. This
% correlation refers to the ROI showing the "largest brain-wide
% association". One can change this value to simulate also the other effect
% sizes reported in Figure 1 (e.g., 0.06 or 0.1).
true_r = 0.2;

% Number of ROIs analyzed in the Marek et al., 2022 study
n_rois = 333;

% Overall number of participants
n_subs = 50000;

% Generate random cognitive ability scores each element simulate cognitive
% performance of an individual (as in brain-behavior correlation studies)
cognitive_ability = randn(n_subs,1);

% Generate random resting state data for the largest brain-wide
% association. As cognitive abilities scores are randomly generated from
% a standard normal distribution, to create correlated RSFC we simply
% multiply cognitive ability values by the pearson r and add gaussian
% noise.
resting_state = true_r.*cognitive_ability + randn(n_subs,1);

% An array storing the number of observations in each simulated sample
% size. These are the same values used in Figure 1
n_obs = [25, 33, 50, 70, 100, 135, 200, 265,...
    375, 525, 725, 1000, 1430, 2000, 2800, 3604];

% Number of resampling. The same value used in the paper
n_boostrap = 1000;

% Significance level. This is the classic *uncorrected* alpha threshold
alpha_level = 0.05;

% Preallocate an array to store simulated brain-behavior correlations at
% different sample sizes and for each random resampling
simulated_corr = nan(n_boostrap,numel(n_obs));

% Preallocate an array to store the level of significane for the 
% brain-behavior correlations at different sample sizes and for each random
% resampling
simulated_p = nan(n_boostrap,numel(n_obs));

% For each random resampling
for b = 1:n_boostrap
    
    % Randomize the order of simulated participants
    shuffle_array = randperm(n_subs);
    
    % For each sample size
    for o = 1:numel(n_obs)
        
        % The selected sample size
        obs = n_obs(o);
        
        % Randomly select observations based on sample size. Please note we
        % are adding to the smallest group of participants (n = 25) more 
        % and more individuals for each bootstrap based on the selected
        % number of observations. This is because we don't want to 
        % introduce a random bias across sample sizes pertaining to the
        % same resampling
        s = shuffle_array(1:obs);
        
        % Select random simulated resting state data
        resting_state_s = resting_state(s);
        
        % Select random simulated cognitive ability data
        cognitive_ability_s = cognitive_ability(s);
        
        % Estimate the correlation and its significance and store them into
        % separate matrices. Please note, here we computed statistical
        % significance based on parametric testing because of two reasons:
        % (1) data are generated from normal distributions, thus
        % assumptions for parametric testing are very likely to be met; (2)
        % if the assumptions are not met (because of any reason), we want
        % to be prone to false positives, as we are interested in getting
        % the number of small sample size studies (i.e., resamplings)
        % reporting an opposite effect as compared to the true one.
        [simulated_corr(b,o),simulated_p(b,o)] = ...
            corr(resting_state_s,cognitive_ability_s);                
        
        
    end
    
end

% Compute the average effect across random resamplings and sample sizes
avg_effect = mean(simulated_corr,1);

% Compute the 2.5th percentile of the effect across random resamplings and 
% sample sizes
lower_prctile_95_effect = prctile(simulated_corr,2.5,1);

% Compute the 97.5th percentile of the effect across random resamplings and 
% sample sizes
upper_prctile_95_effect = prctile(simulated_corr,97.5,1);

% Compute the 0.5th percentile of the effect across random resamplings and 
% sample sizes
lower_prctile_99_effect = prctile(simulated_corr,0.5,1);

% Compute the 99.5th percentile of the effect across random resamplings and 
% sample sizes
upper_prctile_99_effect = prctile(simulated_corr,99.5,1);

% Compute the max of the effect across random resamplings and sample sizes
max_effect = max(simulated_corr,[],1);

% Compute the min of the effect across random resamplings and sample sizes
min_effect = min(simulated_corr,[],1);

% x values for plotting max-min as a shaded area
shaded_x=[n_obs,fliplr(n_obs)];

% y values for plotting max-min as a shaded area
shaded_y=[max_effect,fliplr(min_effect)];

% Preallocate a matrix to store the number of significant positive (i.e.,
% first row) and negative (i.e., second row) associations according to the
% alpha level. Columns are the 16 different sample sizes used in Marek et
% al., 2022 paper. Please note the true effect (i.e., true_r) is positive, 
% thus we should not observe any significant pvalue for negative (spurious)
% correlations. The number of positive/negative significant relationships
% should be interpreted bearing in mind the number of random resamplings
% (i.e., 1000).
sign_p = nan(2,numel(n_obs));

% As sign_p but for storing pvalues corrected for multiple comparisons
sign_pcorr = nan(2,numel(n_obs));

% An array to store the number of negative (i.e., sign-flipped)
% correlations.
n_neg_corr_all = nan(1,numel(n_obs));

% For each sample size
for o = 1:numel(n_obs)
   
    % Identify positive correlations
    pos_corr = simulated_corr(:,o)>0;
    
    % Identify negative correlations
    neg_corr = simulated_corr(:,o)<0;
    
    % Count how many positive correlations
    n_pos_corr = sum(pos_corr);
    
    % Count how many negative correlations
    n_neg_corr = sum(neg_corr);
    
    % Store the number of negative correlations
    n_neg_corr_all(o) = n_neg_corr;
    
    % Select pvalues passing the significance threshold for positive
    % correlations
    pos_corr_sign_p = simulated_p(pos_corr,o) < alpha_level;
    
    % Select pvalues passing the significance threshold for positive
    % correlations after correcting for multiple comparisons
    pos_corr_sign_pcorr = simulated_p(pos_corr,o) < (alpha_level / n_rois);
    
    % Select pvalues passing the significance threshold for negative
    % correlations
    neg_corr_sign_p = simulated_p(neg_corr,o) < alpha_level;
    
    % Select pvalues passing the significance threshold for negative
    % correlations after correcting for multiple comparisons
    neg_corr_sign_pcorr = simulated_p(neg_corr,o) < (alpha_level / n_rois);    
    
    % Number of significant pvalues for positive correlations
    n_pos_corr_sign_p = sum(pos_corr_sign_p);
    
    % Number of significant pvalues for positive correlations after
    % correcting for multiple comparisons
    n_pos_corr_sign_pcorr = sum(pos_corr_sign_pcorr);
    
    % Number of significant pvalues for negative correlations
    n_neg_corr_sign_p = sum(neg_corr_sign_p);
    
    % Number of significant pvalues for negative correlations after
    % correcting for multiple comparisons
    n_neg_corr_sign_pcorr = sum(neg_corr_sign_pcorr);
    
    % Store the number of significant pvalues for positive correlations   
    sign_p(1,o) = n_pos_corr_sign_p;
    
    % Store the number of significant pvalues for negative correlations
    sign_p(2,o) = n_neg_corr_sign_p;
    
    % Store the number of significant pvalues for positive correlations
    % after adjusting for multiple comparisons
    sign_pcorr(1,o) = n_pos_corr_sign_pcorr;
    
    % Store the number of significant pvalues for negative correlations
    % after adjusting for multiple comparisons    
    sign_pcorr(2,o) = n_neg_corr_sign_pcorr;
                
end

% Create a table of the number of uncorrected p for positive and negative
% correlations
sign_p_table = array2table(sign_p,...
    'VariableNames',strsplit(num2str(n_obs)),...
    'RowNames',{sprintf('Positive r, p<%.2f (n)',alpha_level),...
    sprintf('Negative r, p<%.2f (n)',alpha_level)});

disp(sign_p_table)

% Create a table of the number of bonferroni corrected p for positive and 
% negative correlations
sign_pcorr_table = array2table(sign_pcorr,...
    'VariableNames',strsplit(num2str(n_obs)),...
    'RowNames',{sprintf('Positive r, pBonf<%.2f (n)',alpha_level),...
    sprintf('Negative r, pBonf<%.2f (n)',alpha_level)});

disp(sign_pcorr_table)

%% Figures

% Replicating Marek et al., 2022 Figure 1F with simulated data
figure;
S = fill(shaded_x,shaded_y,...
    [.7 .9 .7],...
    'FaceAlpha',0.3,...
    'EdgeColor','none');  
hold on
A = plot(n_obs, avg_effect,...
    'LineStyle','-',...
    'LineWidth',1.2,...
    'Color',[.2 .8 .2]);
L95 = plot(n_obs, lower_prctile_95_effect,...
    'LineStyle',':',...
    'LineWidth',1,...
    'Color',[.2 .2 .2]);
U95 = plot(n_obs, upper_prctile_95_effect,...
    'LineStyle',':',...
    'LineWidth',1,...
    'Color',[.2 .2 .2]);
L99 = plot(n_obs, lower_prctile_99_effect,...
    'LineStyle','--',...
    'LineWidth',1,...
    'Color',[.2 .2 .2]);
U99 = plot(n_obs, upper_prctile_99_effect,...
    'LineStyle','--',...
    'LineWidth',1,...
    'Color',[.2 .2 .2]);
set(gca, 'XScale', 'log',...
    'XTick',[25 200 2000],...
    'XTickLabel',{'25','200','2000'},...
    'YTick',-1:0.5:1,...
    'TickDir','out',...
    'FontName','Arial',...
    'FontSize',18,...
    'YGrid','on');
xlim([min(n_obs),max(n_obs)])
ylim([-1 1])
box off
set(gcf,'color',[1 1 1])
ylabel('Correlation (r)')
xlabel('Sample size')
legend([S,A,L95,L99],...
    {'Min-Max','Avg','95CI','99CI'},'location','southeast','Box','off')
title('Marek et al., 2022 - Figure 1F','FontWeight','normal',...
    'FontName','Arial','FontSize',20)


% Plot number of significant pvalues for correlations of expected (i.e., 
% positive) and opposite (i.e., negative) sign. Results not corrected and
% corrected (Bonferroni method) for multiple comparisons.
figure;
PUNC = plot(n_obs,sign_p(1,:)./n_boostrap,...
    'Color',[0.9843,0.5020,0.4471],...
    'LineWidth',2);
hold on
NUNC = plot(n_obs,sign_p(2,:)./n_boostrap,...
    'Color',[0.5020,0.6941,0.8275],...
    'LineWidth',2);

PCOR = plot(n_obs,sign_pcorr(1,:)./n_boostrap,...
    'Color',[0.9843,0.5020,0.4471],...
    'LineWidth',2,...
    'LineStyle','--');

NCOR = plot(n_obs,sign_pcorr(2,:)./n_boostrap,...
    'Color',[0.5020,0.6941,0.8275],...
    'LineWidth',2,...
    'LineStyle','--');

xlim([min(n_obs),max(n_obs)])
ylim([0 1])
box off
set(gca, 'XScale', 'log',...    
    'XTick',[25 200 2000],...
    'XTickLabel',{'25','200','2000'},...
    'YTick',0:0.2:1,...
    'TickDir','out',...
    'FontName','Arial',...
    'FontSize',18,...
    'YGrid','on');
set(gcf,'color',[1 1 1])
ylabel({'Proportion of','Significant pvalues'})
xlabel('Sample size')

legend([PUNC,NUNC,PCOR,NCOR],...
    {sprintf('p<%.2f | r>0',alpha_level),...
    sprintf('p<%.2f | r<0',alpha_level),...
    sprintf('p_B_o_n_f<%.2f | r>0',alpha_level),...
    sprintf('p_B_o_n_f<%.2f | r<0',alpha_level)},...
    'location','eastoutside','Box','off',...
    'Orientation','Vertical',...
    'FontSize',14)
