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
figure_figure_font_size = 12;

% Set the number of individuals in the population
n_population = 100000;

% Set the magnitude of the effect size in the population
rhos_population = [0.05 0.1 0.15 0.2];

% The minimum and maximum sample size in smaller n studies
sample_size_range = [25 25000];

% Number of log-spaced sample sizes to be tested
n_samples = 32;

% Number of alpha levels to simulate different degrees of correction for
% multiple comparisons
alpha_levels = [5e-2 5e-3 5e-4 5e-5];

% Number of smaller n studies to be simulated for each sample size
n_studies = 1000;

% Number of resamplings to estimate the confidence intervals of the effect
% size (i.e., correlation) in each smaller n study
n_boots = 1000;

% The confidence intervals to be computed
ci_percentiles = 0:0.5:100;

%% Compute sample size adjustment curves

% Number of population effect sizes to be tested
n_rhos = numel(rhos_population);

% Number of alpha levels to be tested
n_alphas = numel(alpha_levels);

% Compute sample size adjustment curves as a function of the effect size in
% the population, the sample size in smaller n study and the level of
% correction for multiple comparisons
sample_size_adjustment = ...
    compute_sample_size_adjustment(rhos_population, ...
    sample_size_range, ...
    n_samples, ...
    alpha_levels, ...
    n_population, ...
    n_studies, ...
    n_boots, ...
    ci_percentiles, ...
    random_seed);

% Plot distribution of actual and simulated correlation values
if startsWith(lower(save_workspace),'y')
    
    % Create workspace filename including path to output directory
    workspace_filename = strcat(output_folder,...
        '/Sample_Size_Adjustment.mat');
    
    % Save workspace
    save(workspace_filename,'-v7.3')
    
end

%% Display figures

% Plot sample size adjustment figures
if startsWith(lower(display_figures),'y')
    
    % Create a colormap using color brewer
    rho_colors = flipud(cbrewer('qual','Set1',n_rhos+1));
    
    % Modify the colormap
    rho_colors = flipud(rho_colors(1:n_rhos,:));
    
    % Define the sample size ticks
    sample_size_ticks = [25 250 2500 25000];
    
    % Define the size of figure in pixels
    figure_size = [680 537 600 600];
    
    % Size of the marker for the effect size in the population
    marker_size = 200;
    
    % Thickness of the line used to plot the gamma distribution
    line_width = 2;
    
    % Size of font in plot
    figure_font_size = 18;
    
    % Set transparency of the histogram
    hist_face_alpha = 0.4;
    
    % Select the id of the alpha value. Plots refer to this alpha level
    id_selected_alpha = 1;
    
    % Get the value of the selected alpha
    selected_alpha = alpha_levels(id_selected_alpha);
    
    % Range for the estimation of the gamma distribution. We explore the
    % range between 0 and 1 as in the simulation correlations are positive
    gamma_xi = 0:0.001:1;
    
    % Determine which sample sizes to plot
    id_selected_samples = [15 20];
    
    % Number of sample sizes to plot
    n_samples_to_plot = numel(id_selected_samples);
    
    % Create a new figure
    figure_title = sprintf('p < %.5f',selected_alpha);
    figure('Position', figure_size, 'Name',figure_title);
    axes('Position',[.25 .25 .45 .45])
    
    hold on
    
    % For each correlation in the population
    for r = 1:n_rhos
        
        % Plot the position of the population effect size in the
        % confidence interval of smaller n studies. As the effect size
        % reported in smaller n studies is inflated (because of the
        % significance threshold), we expect the true effect size to be
        % closer to the lower boundary of the confidence interval for
        % very small n studies and to approach 50th percentile with
        % larger n
        plot(sample_size_adjustment.sample_sizes,...
            sample_size_adjustment.mean_optimal_ci_percentile(r,...
            :,...
            id_selected_alpha),...
            'LineWidth',line_width,'Color',rho_colors(r,:))
        
    end
    
    % Modify the plot appearence
    set(gca, 'XScale', 'log',...
        'XTick',sample_size_ticks,...
        'YTick',0:10:50,...
        'TickDir','out',...
        'FontName','Arial',...
        'FontSize',figure_font_size,...
        'YGrid','on');
    xlim([min(sample_size_adjustment.sample_sizes),...
        max(sample_size_adjustment.sample_sizes)])
    ylim([0 55])
    box off
    set(gcf,'color',[1 1 1])
    ylabel('CI Percentile (%)')
    xlabel('Sample size')
    legend(strcat('r = ',{' '},strsplit(num2str(rhos_population),' ')),...
        'location','southeast','Box','off')
    
    % Save figure
    if startsWith(lower(save_figures),'y')
        
        image_filename = strcat(output_folder,'/',...
            sprintf('Sample_Size_confidence_interval_alpha_%.5f.png',...
            selected_alpha));
        export_fig(image_filename,'-m10','-nocrop', '-transparent','-silent')
        
    end
    
    % For each sample size
    for s = 1:n_samples_to_plot
        
        % Select the sample size
        selected_sample = id_selected_samples(s);
        
        % Create a new figure
        figure_title = sprintf('sample size: %d, p < %.5f',...
            sample_size_adjustment.sample_sizes(selected_sample),...
            selected_alpha);
        figure('Position', figure_size, 'Name',figure_title);
        axes('Position',[.25 .25 .45 .45])
        hold on
        
        % Preallocate graph object
        P = gobjects(n_rhos,1);
        
        % For each effect size in the population
        for r = 1:n_rhos
            
            % Get correlation values for the selected effect size, sample size
            % and significance level
            corr_vals = ...
                squeeze(sample_size_adjustment.all_corr_single_study(r,...
                selected_sample,...
                sample_size_adjustment.all_p_single_study(r,...
                selected_sample,:) < selected_alpha));
            
            % Remove sign flips
            corr_vals = corr_vals(sign(corr_vals) ...
                == sign(rhos_population(r)));
            
            % Fit gamma distribution to the correlation values
            corr_vals_gamma = fitdist(corr_vals,'gamma');
            
            % Get gamma parameters
            gamma_a = corr_vals_gamma.a;
            gamma_b = corr_vals_gamma.b;
            
            % Estimate gamma pdf in the xi range
            corr_vals_gamma_pdf = gampdf(gamma_xi,gamma_a,gamma_b);
            
            % Plot the histogram of correlation values
            histogram(corr_vals,...
                'normalization','pdf',...
                'facecolor',rho_colors(r,:),...
                'FaceAlpha',hist_face_alpha,...
                'edgecolor','none',...
                'binmethod','fd')
            
            % Plot the gamma probability density function
            P(r) = plot(gamma_xi,...
                corr_vals_gamma_pdf,...
                'LineWidth',line_width,...
                'Color',rho_colors(r,:));
            
        end
        
        % For each correlation
        for r = 1:n_rhos
            
            % Plot the effect size in the population
            scatter(rhos_population(r),...
                0,...
                marker_size,...
                'Marker','^',...
                'MarkerFaceColor',rho_colors(r,:),...
                'MarkerEdgeColor',[.2 .2 .2])
            
            % Change plot appearence
            set(gca,...
                'TickDir','out',...
                'FontName','Arial',...
                'FontSize',figure_font_size,...
                'YGrid','on');
            xlim([0, 0.5]);
            box off
            set(gcf,'color',[1 1 1])
            ylabel('Probability density')
            xlabel('Correlation (r)')
            
        end
        
        % Add legend
        legend(P,...
            strcat('r = ',{' '},...
            strsplit(num2str(rhos_population),' '),' '),...
            'location','northeast','Box','off')
        
        % Save figure
        if startsWith(lower(save_figures),'y')
            
            image_filename = strcat(output_folder,'/',...
                sprintf('Sample_Size_effect_inflation_n_%d_alpha_%.5f.png',...
                sample_size_adjustment.sample_sizes(selected_sample),...
                selected_alpha));
            export_fig(image_filename,'-m10','-nocrop', '-transparent','-silent')
            
        end
        
    end
    
end






