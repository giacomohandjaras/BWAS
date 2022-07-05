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

%% Variables

% The four datasets to be imported
path_to_data_rsfc_components = 'raw_data/Component_rsfMRI_Cognition.csv';
path_to_data_rsfc_networks = 'raw_data/Network_rsfMRI_Cognition.csv';
path_to_data_thck_rois = 'raw_data/ROI_Thickness_Cognition.csv';
path_to_data_thck_networks = 'raw_data/Network_Thickness_Cognition.csv';

% The full sample in resting state data
full_sample_size_rsfc = 3604;

% The full sample in cortical thickness
full_sample_size_thck = 3928;

% Precision of the correlation
tolerance = 1e-4;

% Alpha level
alpha = 0.05;

% Range of sample size
sample_size = [25 3000];

% Number of sample sizes at which resamplings are computed
n_samples = 100;

% Number of resamplings
n_boots = 1000;

% Reference value for statistical power
reference_power = 80;

% Refenrence value for sample size. Is just one correlation reaching 80%
% power with less than 1000 individuals?
reference_sample = 1000;

% Seed for generating random numbers
fixed_seed = 15012018;

% Figure size
figure_size = [680 537 600 600];

% Figure font size
figure_font_size = 18;

% Range of the xaxis
xaxis_range = [0 2000];

% X axis ticks
xticks_plot = min(xaxis_range):500:max(xaxis_range);

% Color of the line for each of the four scenarios tested: rsFC components,
% rsFC networks, thickness ROIs, thickness network
plot_colors = [138,201,38;...
255, 89, 94;...
25, 130, 196;...
255, 202, 58]./255;

% Save figure?
save_figures = 'no';

% Where the figure will be stored
output_folder = 'derivatives';

% Fix the seed for reproducibility
rng(fixed_seed);

%% Generate random correlated data

% rsFC components
simulated_data_rsfc_components = ...
    generate_simulated_data(path_to_data_rsfc_components, ...
    full_sample_size_rsfc, ...
    tolerance);

% rsFC networks
simulated_data_rsfc_networks = ...
    generate_simulated_data(path_to_data_rsfc_networks, ...
    full_sample_size_rsfc, ...
    tolerance);

% Thickness ROIs
simulated_data_thck_rois = ...
    generate_simulated_data(path_to_data_thck_rois, ...
    full_sample_size_thck, ...
    tolerance);

% Thickness networks
simulated_data_thck_networks = ...
    generate_simulated_data(path_to_data_thck_networks, ...
    full_sample_size_thck, ...
    tolerance);

%% Resamplings

% rsFC components
bootstrap_output_rsfc_components = ...
    bootstrap_correlation(simulated_data_rsfc_components.behavior, ...
    simulated_data_rsfc_components.brain, ...
    alpha, ...
    sample_size, ...
    n_samples, ...
    n_boots);

% rsFC networks
bootstrap_output_rsfc_networks = ...
    bootstrap_correlation(simulated_data_rsfc_networks.behavior, ...
    simulated_data_rsfc_networks.brain, ...
    alpha, ...
    sample_size, ...
    n_samples, ...
    n_boots);

% Thickness ROIs
bootstrap_output_thck_rois = ...
    bootstrap_correlation(simulated_data_thck_rois.behavior, ...
    simulated_data_thck_rois.brain, ...
    alpha, ...
    sample_size, ...
    n_samples, ...
    n_boots);

% Thickness networks
bootstrap_output_thck_networks = ...
    bootstrap_correlation(simulated_data_thck_networks.behavior, ...
    simulated_data_thck_networks.brain, ...
    alpha, ...
    sample_size, ...
    n_samples, ...
    n_boots);

%% Estimate power

% rsFC components
statistical_power_rsfc_components = ...
    compute_power(bootstrap_output_rsfc_components);

% rsFC networks
statistical_power_rsfc_networks = ...
    compute_power(bootstrap_output_rsfc_networks);

% Thickness ROIs
statistical_power_thck_rois = ...
    compute_power(bootstrap_output_thck_rois);

% Thickness networks
statistical_power_thck_networks = ...
    compute_power(bootstrap_output_thck_networks);

%% Count FDR corrected powered correlations

% rsFC components 
rsfc_components_n_powered_corr = ...
    sum(statistical_power_rsfc_components.fdr >= reference_power,2);

% rsFC networks
rsfc_networks_n_powered_corr = ...
    sum(statistical_power_rsfc_networks.fdr >= reference_power,2);

% Thickness ROIs
thck_rois_n_powered_corr = ...
    sum(statistical_power_thck_rois.fdr >= reference_power,2);

% Thickness networks
thck_networks_n_powered_corr = ...
    sum(statistical_power_thck_networks.fdr >= reference_power,2);

%% Plot results
    
figure('Position', figure_size, 'Name','FDR corrected');
axes('Position',[.25 .25 .45 .45])

hold on

plot(bootstrap_output_rsfc_components.sample_sizes,...
    rsfc_components_n_powered_corr,...
    'LineWidth',2.5,...
    'Color',plot_colors(1,:))

plot(bootstrap_output_rsfc_networks.sample_sizes,...
    rsfc_networks_n_powered_corr,...
    'LineWidth',2.5,...
    'Color',plot_colors(2,:))

plot(bootstrap_output_thck_rois.sample_sizes,...
   thck_rois_n_powered_corr,...
   'LineWidth',2.5,...
   'Color',plot_colors(3,:))

plot(bootstrap_output_thck_networks.sample_sizes,...
   thck_networks_n_powered_corr,...
   'LineWidth',2.5,...
   'Color',plot_colors(4,:))

set(gca,'XTick',xticks_plot,...
    'YTick',[4 8 12 16 20 24],...
    'TickDir','out',...
    'FontName','Arial',...
    'FontSize',figure_font_size,...
    'YGrid','on',...
    'XGrid','on');

xlim(xaxis_range)
ylabel({'FDRc Correlations','≥ 80% power (n)'})
xlabel('Sample size (n)')
axis square
box off
set(gcf,'color',[1 1 1])
legend({'rsFC components',...
    'rsFC networks',...
    'thck ROIs',...
    'thck networks'},...
    'Location','northwest',...
    'box','off')

if startsWith(lower(save_figures),'y')
    
    image_filename = strcat(output_folder,'/count_powered_correlations.png');
    export_fig(image_filename,'-m10','-nocrop', '-transparent','-silent')
    
end



