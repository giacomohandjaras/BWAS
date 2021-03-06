function plot_correlation_variability(bootstrap_output, selected_percentiles, figure_size, plot_colormap, sample_size_ticks, figure_font_size, xaxis_type, dataset_filename, save_figures, output_folder)
% coded by Luca Cecchetti & Giacomo Handjaras, IMT-Lucca, Italy
% vers 20220619

%figure_title = strrep(strrep(dataset_filename,'.csv',''),'_',' ');

% The best correlation in the full sample
[~, id_max_abs_corr] = max(abs(bootstrap_output.r_full));

% Select data for the best correlation in the full sample
best_corr_bootstrap = ...
    bootstrap_output.simulated_corr(:,:,id_max_abs_corr);

% Nth percentile across bootstraps for each sample size
best_corr_bootstrap_percentiles = ...
    prctile(best_corr_bootstrap,selected_percentiles);

% Average correlation across bootstraps for each sample size
best_corr_bootstrap_average = mean(best_corr_bootstrap);

% x values for plotting max-min as a shaded area
shaded_x=[bootstrap_output.sample_sizes,...
    fliplr(bootstrap_output.sample_sizes)];

% y values for plotting max-min as a shaded area
shaded_y=[best_corr_bootstrap_percentiles(end,:),...
    fliplr(best_corr_bootstrap_percentiles(1,:))];

figure('Position', figure_size, 'Name','Sampling variability');
axes('Position',[.25 .25 .45 .45])
line([min(bootstrap_output.sample_sizes),...
    max(bootstrap_output.sample_sizes)],...
    [0,0],'LineWidth',1.2,...
    'LineStyle','--',...
    'Color',[.9 .2 .2]);
hold on

S = fill(shaded_x,shaded_y,...
    plot_colormap(1,:),...
    'FaceAlpha',0.3,...
    'EdgeColor','none');

A = plot(bootstrap_output.sample_sizes, best_corr_bootstrap_average,...
    'LineStyle','-',...
    'LineWidth',1.5,...
    'Color',plot_colormap(1,:));

P1 = plot(repmat(bootstrap_output.sample_sizes,2,1)',...
    best_corr_bootstrap_percentiles([2,end-1],:)',...
    'LineStyle','--',...
    'LineWidth',1,...
    'Color',[.2 .2 .2]);

P2 = plot(repmat(bootstrap_output.sample_sizes,2,1)',...
    best_corr_bootstrap_percentiles([3,end-2],:)',...
    'LineStyle',':',...
    'LineWidth',1,...
    'Color',[.2 .2 .2]);

set(gca, 'XScale', xaxis_type,...
    'XTick',sample_size_ticks,...
    'XTickLabel',strsplit(num2str(sample_size_ticks),' '),...
    'YTick',-1:0.5:1,...
    'TickDir','out',...
    'FontName','Arial',...
    'FontSize',figure_font_size,...
    'YGrid','on');

xlim([min(bootstrap_output.sample_sizes),...
    max(bootstrap_output.sample_sizes)])
ylim([-1 1])

box off
axis square
set(gcf,'color',[1 1 1])
ylabel('Correlation (r)')
xlabel('Sample size')
legend([S,A,P1(1),P2(1)],...
    {'Min-Max','??','99CI','95CI'},'location','southeast','Box','off')
pbaspect([1 1 1])


% title(figure_title,'FontWeight','normal',...
%     'FontName','Arial','FontSize',figure_font_size)

if startsWith(lower(save_figures),'y')
    
    image_filename = strcat(output_folder,'/',...
        strrep(dataset_filename,'.csv','_sampling_variability.png'));
    export_fig(image_filename,'-m10','-nocrop', '-transparent','-silent')
    
end

end
