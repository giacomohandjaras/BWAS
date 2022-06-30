%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% coded by Luca Cecchetti & Giacomo Handjaras, IMT-Lucca, Italy
%%% vers 20220630
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
									

% Import data. CSV are built starting from figure 4 source data of M&TC
% Nature paper. The original figure depicts results for the association
% between the first canonical variate obtained from rsFC (or thickness) and
% the one obtained from NIH Toolbox (or CBCL).
% We focus on replication values (i.e., out-of-sample correlations) because
% they are less inflated than discovery correlations.

% Determine whether figures should be displayed on screen
display_figures = 'yes';

% Determine whether figures should be saved on disk
save_figures = 'no';

% Where the figure will be stored
output_folder = 'derivatives';

% Import CCA results for rsFC ~ cognition
cca_rsfmri_cognition = ...
    csvread('raw_data/CCA_rsfMRI_Cognition_Replication.csv');

% Import CCA results for rsFC ~ psychopathology
cca_rsfmri_psychopathology = ...
    csvread('raw_data/CCA_rsfMRI_Psychopathology_Replication.csv');

% Import CCA results for cortical thickness ~ cognition
cca_thickness_cognition = ...
    csvread('raw_data/CCA_Thickness_Cognition_Replication.csv');

% Import CCA results for cortical thickness ~ psychopathology
cca_thickness_psychopathology = ...
    csvread('raw_data/CCA_Thickness_Psychopathology_Replication.csv');

% Catenate datasets on the third dimension
cca_all_data = cat(3,cca_rsfmri_cognition,...
    cca_rsfmri_psychopathology,...
    cca_thickness_cognition,...
    cca_thickness_psychopathology);

% Samples at which M&TC evaluated the performance of CCA. Numbers are
% obtained from the abcd_cca.m script found here:
% https://gitlab.com/DosenbachGreene/bwas/-/blob/main/abcd_cca.m
% We omit the last sample size as it depends on the brain imaging
% technique. The remaining sample sizes are identical across the
% four conditions explored.
samples = [25, 33, 45, 60, 80, 100, 145, 200,...
    265, 350, 460, 615, 825, 1100, 1475]';

% Reference power
reference_power = 80;

% Level of significance for the computation of statistical power
alpha = 0.05;

% Determine the number of samples
n_samples = numel(samples);

% Determine the number of resamplings at each sample size. For the
% multivariate analysis M&TC performed 100 resamplings
n_resamplings = size(cca_all_data,2);

% Number of conditions
n_conditions = size(cca_all_data,3);

% Preallocate matrix to store pvalues
pvalues = nan(n_samples,n_resamplings,n_conditions);

% For each condition
for c = 1:n_conditions
    
    % Feedback about progress
    fprintf('Processing dataset %d of %d\n',c,n_conditions);
    
    % Extract correlations obtained at all resamplings and sample sizes. We
    % take the absolute value as we are interested to significance only
    cancorr_r = abs(cca_all_data(:,:,c));
    
    % For each sample size
    for i = 1:n_samples
        
        % For each resampling
        for j = 1:n_resamplings
            
            % Extract correlation value
            r = cancorr_r(i,j);
            
            % Transform correlation to t-value
            t = r*sqrt((samples(i)-2)/(1-r^2));
            
            % Estimate statistical significance. We test negative tail
            % because conversion to pvalue yields more accurate values.
            % Otherwise very strong effects will get pvalue = 0
            pvalues(i,j,c) = tcdf(-t,samples(i))*2;
            
        end
        
    end
    
end

% Compute statistical power.
statistical_power = ...
    squeeze(...
    sum(pvalues < alpha,2))...
    ./n_resamplings.*100;

% Convert matrix to a table to improved readability
statistical_power_table = ...
    array2table(cat(2,samples,statistical_power),...
    'VariableNames',{'Sample size',...
    'power: rsFC~cognition',...
    'power: rsFC~psych',...
    'power: thick~cognition',...
    'power: thick~psych'});

%% Plot CCA power

% Plot multivariate power
if startsWith(lower(display_figures),'y')
    
    % Define plot variables
    sample_size_ticks = [25 50 100 250 550 1250];
    figure_size = [680 537 600 600];
    figure_font_size = 18;
    xaxis_type = 'log';
    xmin = 15;
    line_color = [51,160,44;...
        106,61,154;...
        178,223,138;...
        202,178,214]./255;
    
    % Plot
    figure('Position', figure_size, 'Name','CCA');
    axes('Position',[.25 .25 .45 .45]);
    hold on
    
    for c = 1:n_conditions
        
        plot(samples,statistical_power(:,c),...
            'LineStyle','-',...
            'LineWidth',2,...
            'Color',line_color(c,:));
        
    end
    
    line([xmin,...
        max(samples)],...
        [reference_power,...
        reference_power],...
        'LineStyle','--',...
        'LineWidth',1,...
        'Color',[.2 .2 .2]);
    
    legend({'rsFC~Cog',...
        'rsFC~Psy',...
        'thck~Cog',...
        'thck~Psy'},...
        'Location','west',...
        'box','off')
    
    set(gca, 'XScale', xaxis_type,...
        'XTick',sample_size_ticks,...
        'YTick',0:20:100,...
        'TickDir','out',...
        'FontName','Arial',...
        'FontSize',figure_font_size,...
        'YGrid','on',...
        'XGrid','on');
    xlim([xmin,...
        max(samples)])
    ylim([0 100])
    box off
    set(gcf,'color',[1 1 1])
    ylabel('Power (%)')
    xlabel('Sample size')
    
    if startsWith(lower(save_figures),'y')
        
        image_filename = strcat(output_folder,'/cca_replication_power.png');
        export_fig(image_filename,'-m10','-nocrop', '-transparent','-silent')
        
    end            
    
end


