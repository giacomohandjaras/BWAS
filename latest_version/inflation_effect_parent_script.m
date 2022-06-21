%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% coded by Luca Cecchetti & Giacomo Handjaras, IMT-Lucca, Italy
%%% vers 20220619
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc

%% External libraries & code

% Our functions
addpath('code/')

% Extra libraries & functions 
addpath('external_code/')   % fdr_bh.m: 
							% Groppe, D.M., Urbach, T.P., & Kutas, M. (2011) 
							% Mass univariate analysis of event-related 
							% brain potentials/fields I: A critical tutorial
							% review. Psychophysiology, 48(12) pp. 1711-1725,
							% DOI: 10.1111/j.1469-8986.2011.01273.x
addpath('external_code/export_fig/')    % export_fig.m: 
										% Copyright (c) 2014, 
										% Oliver J. Woodford, 
										% Yair M. Altman
addpath('external_code/cbrewer/')   % cbrewer.m
									% Copyright (c) 2015, Charles Robert

%% Define simulation parameters

% Fix the seed for reproducibility
random_seed = 15012018;

% Folder to store analysis output
output_folder = 'derivatives';

% Determine whether results should be written to disk
save_workspace = 'no';

% Determine whether figures should be displayed
display_figures = 'yes';

% Determine whether figures should be written to disk
save_figures = 'no';

% Size of font in figures
figure_font_size = 12;

% Set the number of individuals in the full sample
full_sample_size = 4000;

% Set the magnitude of the effect size in the full sample
rhos_population = 0.05:0.05:0.2;

n_rhos = numel(rhos_population);

% Tolerance in the absolute difference between original correlation values
% and simulated correlation values
tolerance = 1e-4;

% Significance level. This threshold is used to estimate significant
% uncorrected, FDR corrected and Bonferroni corrected results
alpha = 0.05;

% The range in sample size for estimating the bootstrapped correlations.
% The minimum sample size is identical to the one used in the original
% Nature publication
sample_size = [25 3000];

% Number of sample sizes used in the computation of bootstrapped
% correlations
n_samples = 100;

% Number of random resamplings
n_boots = 1000;

%% Run analyses - Generate simulated data (Updated)
% Fix the seed for reproducibility
rng(random_seed);

% Generate random correlated data
simulated_data = ...
    generate_simulated_data_single_corr(rhos_population, ...
    full_sample_size, ...
    tolerance);

%% Run analyses - Bootstrapped correlations (Updated)
% Run bootstrap
bootstrap_output = bootstrap_correlation(simulated_data.behavior, ...
    simulated_data.brain, ...
    alpha, ...
    sample_size, ...
    n_samples, ...
    n_boots);

%% Type M error

% Estimate type M errors
type_m_errors = ...
    compute_type_m_alternative(bootstrap_output,'mean');

% Estimate the critical r given a specific combination of significance
% level and sample size.
r_crit_uncorrected = ...
    transform_p_to_r(alpha,1,cat(2,...
    bootstrap_output.sample_sizes,...
    full_sample_size));

% Inflation rate for the critical r when threshold is uncorrected
inflation_percentage_r_crit_uncorrected = ...
    (r_crit_uncorrected(1:end-1) ...
    ./ r_crit_uncorrected(end) - 1) ...
    .*100;

% Add inflation rate of critical uncorrected r to the type_m_errors
% structure
type_m_errors.inflation_percentage_r_crit_uncorrected = ...
    inflation_percentage_r_crit_uncorrected;

% Create a colormap using color brewer
rho_colors = flipud(cbrewer('qual','Set1',n_rhos+1));

% Modify the colormap
rho_colors = flipud(rho_colors(1:n_rhos,:));

% Define the sample size ticks
sample_size_ticks = [25 250 2500 25000];

% Define the size of figure in pixels
figure_size = [680 537 600 600];

% Thickness of the line used to plot the gamma distribution
line_width = 2;

% Create a new figure
figure('Position', figure_size, 'Name','Inflation');
axes('Position',[.25 .25 .45 .45])

hold on

% For each correlation value
for r = 1:n_rhos
    
    % Plot inflation rate as a function of sample size
    plot(bootstrap_output.sample_sizes,...
        type_m_errors.inflation_percentage_thresholded_uncorrected(:,r),...
        'LineWidth',line_width,'Color',rho_colors(r,:))
    
end

% Plot inflation rate for expected correlation at alpha level as a function
% of effect size
plot(bootstrap_output.sample_sizes,...
    type_m_errors.inflation_percentage_r_crit_uncorrected,...
    'LineWidth',line_width,'Color',[.2 .2 .2],...
    'LineStyle',':')

% Define legend labels
legend_labels = cat(2,strcat('r = ',{' '},strsplit(num2str(rhos_population),' ')),...
    sprintf('r_c_r_i_t \\alpha = %.2f',alpha));

% Modify the plot appearence
set(gca, 'XScale', 'log',...
    'XTick',sample_size_ticks,...
    'TickDir','out',...
    'FontName','Arial',...
    'FontSize',figure_font_size,...
    'YGrid','on');
xlim([min(bootstrap_output.sample_sizes),...
    max(bootstrap_output.sample_sizes)])
box off
set(gcf,'color',[1 1 1])
ylabel('Inflation (%)')
xlabel('Sample size')
legend(legend_labels,...
    'location','northeast','Box','off')

% Save figure
if startsWith(lower(save_figures),'y')
    
    image_filename = strcat(output_folder,'/Inflation_effect.png');
    export_fig(image_filename,'-m10','-nocrop', '-transparent','-silent')
    
end

    
% Check if workspace should be written to disk
if startsWith(lower(save_workspace),'y')
        
    % Create workspace filename including path to output directory
    workspace_filename = strcat(output_folder,'/',...
    'Inflation_effect.mat');
        
    % Save workspace
    save(workspace_filename,'-v7.3')
        
end
