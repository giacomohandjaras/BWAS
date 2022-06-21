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

% Define the output file prefix
output_file_prefix = 'HCPdata_ROI_Thickness_PMAT24_A_CR.csv';

% Folder to store analysis output
output_folder = 'derivatives';

% Significance level. This threshold is used to estimate significant
% uncorrected, FDR corrected and Bonferroni corrected results
alpha = 0.05;

% Value for reference statistical power
reference_power = 80;

% The range in sample size for estimating the bootstrapped correlations.
% The minimum sample size is identical to the one used in the original
% Nature publication
sample_size = [25 1000];

% The range in sample size for estimating the replicability of results in
% bootstrapped correlations when doing the half-split
sample_size_replication = [25 500];

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
save_figures = 'no';

% Size of font in figures
figure_font_size = 12;

%% Prepare HCP data

% Load HCP unrestricted data. We are not allowed to directly share HCP data
% but you can apply to get access to unrestricted data, which is sufficient
% to replicate the findings we report in the response letter. Change the
% name of the csv file accordingly
hcp_data = readtable('raw_data/unrestricted_lucacecchetti_3_21_2022_1_1_51.csv');

% Get the index of columns storing Freesurfer cortical thickness values
idx_brain_data = startsWith(hcp_data.Properties.VariableNames,'FS_') & ...
    endsWith(hcp_data.Properties.VariableNames,'Thck');

% Get the ROI name of Freesurfer cortical thickness data
% roiname_brain_data = hcp_data.Properties.VariableNames(...
%     startsWith(hcp_data.Properties.VariableNames,'FS_') & ...
%     endsWith(hcp_data.Properties.VariableNames,'Thck'));

% Get brain data from the HCP table
brain_data = table2array(hcp_data(:,idx_brain_data));

% Get behavioral data from HCP table. We will focus on PMAT24_A_CR as the
% correlation in the full sample is 0.1444
behav_data = hcp_data.PMAT24_A_CR;

% Aggregate behavioral and brain data into a new table so that it will be
% easier to remove missing values
aggregated_table = cat(2,behav_data,brain_data);

% Remove missing values (participants)
aggregated_table_no_missing = rmmissing(aggregated_table);

% Get the number of participants
full_sample_size = size(aggregated_table_no_missing,1);

% Split behavioral data from brain data after missing values are removed
behav_data_no_missing = aggregated_table_no_missing(:,1);

% Split brain data from behavioral data after missing values are removed
brain_data_no_missing = aggregated_table_no_missing(:,2:end);

% Array to split observations into discovery and replication datasets
splitting_array = randperm(full_sample_size);

% Id of observations to be included in the discovery dataset
discovery_id = splitting_array(1:round(full_sample_size/2));

% Id of observations to be included in the replication dataset
replication_id = splitting_array(round(full_sample_size/2)+1:end);

% Compute the association between brain and behavior in the full sample
r_full = corr(behav_data_no_missing,brain_data_no_missing);

% Prepare a data structure that will serve as input for the
% bootstrap_correlation function
actual_data.behavior = behav_data_no_missing;
actual_data.brain = brain_data_no_missing;
actual_data.discovery_behavior = behav_data_no_missing(discovery_id);
actual_data.discovery_brain = brain_data_no_missing(discovery_id,:);
actual_data.replication_behavior = behav_data_no_missing(replication_id);
actual_data.replication_brain = brain_data_no_missing(replication_id,:);
actual_data.original_corr_values = r_full';
actual_data.simulated_corr_values = r_full';

%% Run analyses - Bootstrapped correlations (Updated)
% Run bootstrap
bootstrap_output = bootstrap_correlation(actual_data.behavior, ...
    actual_data.brain, ...
    alpha, ...
    sample_size, ...
    n_samples, ...
    n_boots);

% Plot distribution of actual and simulated correlation values
if startsWith(lower(display_figures),'y')
    
    % Figure parameters
    plot_colormap = [.5 .5 .5];
    selected_percentiles = [0 0.5 2.5 97.5 99.5 100];
    sample_size_ticks = [25 100 500 2000];
    figure_size = [680 537 600 600];
    xaxis_type = 'log';
    
    % Plot correlation variability
    plot_correlation_variability(bootstrap_output,...
        selected_percentiles,...
        figure_size,...
        plot_colormap,...
        sample_size_ticks,...
        figure_font_size,...
        xaxis_type,...
        output_file_prefix,...
        save_figures,...
        output_folder);
    
end

% Check if at least one correlation pass the Bonferroni correction in the
% full sample. Power, type m, type s, replication rate and fpr are computed
% only if this is the case
if sum(bootstrap_output.p_full_significant_bonferroni) > 0
    
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
        sample_size_ticks = [25 250 900];
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
            output_file_prefix,...
            save_figures,...
            output_folder)
        
    end
    
    %% Run analyses - Type M (Updated)
    % Estimate type M errors as in the original paper (with >50%, >100% and
    % >200% inflation bins)
    % type_m_errors = compute_type_m(bootstrap_output);
    
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
    
    % Estimate the critical r given a specific combination of bonferroni
    % corrected significance level and sample size
    r_crit_bonferroni = ...
        transform_p_to_r(alpha,...
        size(actual_data.original_corr_values,1),...
        cat(2,...
        bootstrap_output.sample_sizes,...
        full_sample_size));
    
    % Inflation rate for the critical r when threshold is uncorrected
    inflation_percentage_r_crit_bonferroni = ...
        (r_crit_bonferroni(1:end-1) ...
        ./ r_crit_bonferroni(end) - 1) ...
        .*100;
    
    % Add inflation rate of critical uncorrected r to the type_m_errors
    % structure
    type_m_errors.inflation_percentage_r_crit_uncorrected = ...
        inflation_percentage_r_crit_uncorrected;
    
    % Add inflation rate of critical bonferroni corrected r to the
    % type_m_errors structure
    type_m_errors.inflation_percentage_r_crit_bonferroni = ...
        inflation_percentage_r_crit_bonferroni;
    
    % Plot Type M
    if startsWith(lower(display_figures),'y')
        
        % Figure parameters
        figure_size = [680 537 600 600];
        sample_size_ticks = [25 250 900];
        alpha_color = 0.25;
        xaxis_type = 'log';
        %colormap_correlations = flipud(cbrewer('qual', 'Pastel1', 3));
        colormap_correlations = repmat([.5 .5 .5],3,1);        
        colormap_correlations = cat(2,colormap_correlations,...
            repelem(alpha_color,3,1));
        %colormap_average = flipud(cbrewer('qual', 'Set1', 3));
        colormap_average = repmat([.3 .3 .3],3,1);
        
        % Plot type M error rate as a function of sample size
        plot_type_m(type_m_errors,...
            bootstrap_output,...
            required_sample_sizes,...
            figure_size,...
            sample_size_ticks,...
            figure_font_size,...
            xaxis_type,...
            colormap_correlations,...
            colormap_average,...
            output_file_prefix,...
            save_figures,...
            output_folder);
        
    end
    
    %% Run analyses - Type S (Updated)
    % Estimate Type S errors
    type_s_errors = compute_type_s(bootstrap_output);
    
    % Plot Type S
    if startsWith(lower(display_figures),'y')
        
        % Figure parameters
        figure_size = [680 537 600 600];
        sample_size_ticks = [25 250 900];
        insert_sample_size_ticks = [75 750];        
        alpha_color = 0.25;
        xaxis_type = 'log';
        %colormap_correlations = flipud(cbrewer('qual', 'Pastel1', 3));
        colormap_correlations = repmat([.5 .5 .5],3,1); 
        colormap_correlations = cat(2,colormap_correlations,...
            repelem(alpha_color,3,1));
        %colormap_average = flipud(cbrewer('qual', 'Set1', 3));
        colormap_average = repmat([.3 .3 .3],3,1);
        
        % Plot type S error rate as a function of sample size
        plot_type_s(type_s_errors,...
            bootstrap_output,...
            required_sample_sizes,...
            figure_size,...
            sample_size_ticks,...
            insert_sample_size_ticks,...
            figure_font_size,...
            xaxis_type,...
            colormap_correlations,...
            colormap_average,...
            output_file_prefix,...
            save_figures,...
            output_folder);
        
    end
    
    %% Run analyses - Type I (Updated)
    % Estimate type I error
    false_positive_rate = compute_fpr(bootstrap_output);
    
    % Plot Type I error
    if startsWith(lower(display_figures),'y')
        
        % Figure parameters
        reference_fpr = 5;
        figure_size = [680 537 600 600];
        sample_size_ticks = [25 250 900];
        alpha_color = 0.25;
        xaxis_type = 'log';
        %colormap_correlations = flipud(cbrewer('qual', 'Pastel1', 3));
        colormap_correlations = repmat([.5 .5 .5],3,1); 
        colormap_correlations = cat(2,colormap_correlations,...
            repelem(alpha_color,3,1));
        %colormap_average = flipud(cbrewer('qual', 'Set1', 3));
        colormap_average = repmat([.3 .3 .3],3,1);
        
        % Plot type I error rate as a function of sample size
        plot_fpr(false_positive_rate,...
            bootstrap_output,...
            figure_size,...
            reference_fpr,...
            sample_size_ticks,...
            figure_font_size,...
            xaxis_type,...
            colormap_correlations,...
            colormap_average,...
            output_file_prefix,...
            save_figures,...
            output_folder);
        
    end
    
    %% Run analyses - Probability of Replication (Updated)
    
    % The best correlation in the full discovery sample
    [~, id_max_abs_corr] = ...
        max(abs(bootstrap_output.r_full));    
    
    % Run bootstrap in discovery sample
    bootstrap_output_discovery = ...
        bootstrap_correlation(actual_data.discovery_behavior, ...
        actual_data.discovery_brain, ...
        alpha, ...
        sample_size_replication, ...
        n_samples, ...
        n_boots);
    
    % Run bootstrap in replication sample
    bootstrap_output_replication = ...
        bootstrap_correlation(actual_data.replication_behavior, ...
        actual_data.replication_brain, ...
        alpha, ...
        sample_size_replication, ...
        n_samples, ...
        n_boots);
    
    % Estimate Replication Probability
    replication_rate = ...
        compute_replication(bootstrap_output_discovery,...
        bootstrap_output_replication,...
        id_max_abs_corr);
    
    % Plot Replication Probability
    if startsWith(lower(display_figures),'y')
        
        % Figure parameters
        figure_size = [680 537 600 600];
        sample_size_ticks = [25 100 500];
        alpha_color = 0.25;
        xaxis_type = 'log';
        %colormap_correlations = flipud(cbrewer('qual', 'Pastel1', 3));
        colormap_correlations = repmat([.5 .5 .5],3,1);
        colormap_correlations = cat(2,colormap_correlations,...
            repelem(alpha_color,3,1));
        %colormap_average = flipud(cbrewer('qual', 'Set1', 3));
        colormap_average = repmat([.3 .3 .3],3,1);
        
        % Plot Replication Probability as a function of sample size
        plot_replication_rate(replication_rate,...
            bootstrap_output_discovery,...
            bootstrap_output.sample_sizes,...
            required_sample_sizes,...
            figure_size,...
            sample_size_ticks,...
            figure_font_size,...
            xaxis_type,...
            colormap_correlations,...
            colormap_average,...
            output_file_prefix,...
            save_figures,...
            output_folder)
        
    end
    
    % If no correlations pass the Bonferroni corrected threshold, then exit
    % without computing metrics and print a feedback
else
    
    % Print feedback
    fprintf('\n')
    fprintf(' +++++++++++++++++++++++++++++++++++++++++++ \n')
    fprintf(' + WARNING                                 + \n')
    fprintf(' + No correlations pass the Bonferroni     + \n')
    fprintf(' + corrected threshold in the full sample. + \n')
    fprintf(' + Exit without computing metrics.         + \n')
    fprintf(' +++++++++++++++++++++++++++++++++++++++++++ \n')
    
end

%% Save data to disk (Updated)

% Check if workspace should be written to disk
if startsWith(lower(save_workspace),'y')
    
    % Create workspace filename including path to output directory
    workspace_filename = strcat(output_folder,'/',...
        strrep(output_file_prefix,'.csv','.mat'));
    
    % Save workspace
    save(workspace_filename,'-v7.3')
    
end










