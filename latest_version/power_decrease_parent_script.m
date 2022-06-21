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

%% Input data

% Filename of the csv storing raw data
dataset_filename = 'Simulation_Power.csv';

%% Define simulation parameters

% Folder storing input data, which are csv files obtained from source data
% of Figure 1 in Marek and Tervo-Clemmens et al., 2022
raw_data_folder = 'raw_data';

% Folder to store analysis output
output_folder = 'derivatives';

% Sample sizes in the full sample
full_sample_sizes = [4000 16000 256000];

% Number of full sample sizes to be tested
n_samples_full = numel(full_sample_sizes);

% Tolerance in the absolute difference between original correlation values
% and simulated correlation values
tolerance = 1e-4;

% Significance level. This threshold is used to estimate significant
% uncorrected, FDR corrected and Bonferroni corrected results
alpha = 0.05;

% Value for reference statistical power
reference_power = 80;

% The range in sample size for estimating the bootstrapped correlations.
% The minimum sample size is identical to the one used in the original
% Nature publication
sample_size = [25 3000];

% Number of sample sizes used in the computation of bootstrapped
% correlations
n_samples = 100;

% Number of random resamplings
n_boots = 1000;

% Seed for generating random numbers
fixed_seed = 15012018;

% Determine whether results should be written to disk
save_workspace = 'no';

% Determine whether figures should be displayed
display_figures = 'yes';

% Determine whether figures should be written to disk
save_figures = 'np';

% Size of font in figures
figure_font_size = 12;

%% Run analyses - Generate simulated data (Updated)
% Define the path to dataset file
path_to_data = strcat(raw_data_folder,'/',dataset_filename);

% For each sample size
for s = 1:n_samples_full
    
    dataset_filename_with_sample = ...
        strcat(strrep(dataset_filename,'.csv',''),...
        '_n_',num2str(full_sample_sizes(s)),'.csv');
    
    % Fix the seed for reproducibility
    rng(fixed_seed);
    
    % Generate random correlated data
    simulated_data = ...
        generate_simulated_data(path_to_data, ...
        full_sample_sizes(s), ...
        tolerance);
    
    % Run bootstrap
    bootstrap_output = bootstrap_correlation(simulated_data.behavior, ...
        simulated_data.brain, ...
        alpha, ...
        sample_size, ...
        n_samples, ...
        n_boots);
    
    %% Run analyses - Power (Updated)
    % Estimate power
    statistical_power = compute_power(bootstrap_output);
    
    % Estimate minimum sample size for reference power considering the best
    % correlation and uncorrected threshold
    [~,reference_power_sample_size_uncorrected] = ...
        min(...
        abs(...
        statistical_power.best_uncorrected - reference_power...
        )...
        );
    
    % Save it in a structure
    required_sample_sizes.uncorrected = ...
        reference_power_sample_size_uncorrected;
    
    % Estimate minimum sample size for reference power considering the best
    % correlation and fdr corrected threshold
    [~,reference_power_sample_size_fdr] = ...
        min(...
        abs(...
        statistical_power.best_fdr - reference_power...
        )...
        );
    
    % Save it in a structure
    required_sample_sizes.fdr = ...
        reference_power_sample_size_fdr;
    
    % Estimate minimum sample size for reference power considering the best
    % correlation and bonferroni corrected threshold
    [~,reference_power_sample_size_bonferroni] = ...
        min(...
        abs(...
        statistical_power.best_bonferroni - reference_power...
        )...
        );
    
    % Save it in a structure
    required_sample_sizes.bonferroni = ...
        reference_power_sample_size_bonferroni;
    
    % Plot power
    if startsWith(lower(display_figures),'y')
        
        % Figure parameters
        figure_size = [680 537 600 600];
        sample_size_ticks = [25 100 500 2000];
        alpha_color = 0.25;
        %colormap_correlations = flipud(cbrewer('qual', 'Pastel1', 3));
        colormap_correlations = repmat([.5 .5 .5],3,1);
        colormap_correlations = cat(2,colormap_correlations,...
            repelem(alpha_color,3,1));
        %colormap_average = flipud(cbrewer('qual', 'Set1', 3));
        colormap_average = repmat([.3 .3 .3],3,1);
        xaxis_type = 'log';
        
        % Plot power as a function of sample size
        plot_power(statistical_power,...
            bootstrap_output,...
            reference_power,...
            required_sample_sizes,...
            figure_size,...
            sample_size_ticks,...
            figure_font_size,...
            xaxis_type,...
            colormap_correlations,...
            colormap_average,...
            dataset_filename_with_sample,...
            save_figures,...
            output_folder)
        
    end
    
    %% Save data to disk (Updated)
    
    % Check if workspace should be written to disk
    if startsWith(lower(save_workspace),'y')
        
        % Create workspace filename including path to output directory
        workspace_filename = strcat(output_folder,'/',...
            strrep(dataset_filename_with_sample,'.csv','.mat'));
        
        % Save workspace
        save(workspace_filename,'-v7.3')
        
    end
    
    
end
