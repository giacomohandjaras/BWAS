function plot_type_s(type_s_errors, bootstrap_output, required_sample_sizes, figure_size, sample_size_ticks, insert_sample_size_ticks, figure_font_size, xaxis_type, colormap_correlations, colormap_average, dataset_filename, save_figures, output_folder)
% coded by Luca Cecchetti & Giacomo Handjaras, IMT-Lucca, Italy
% vers 20220619

%figure_title = strrep(strrep(dataset_filename,'.csv',''),'_',' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UNCORRECTED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position', figure_size, 'Name','Uncorrected');
axes('Position',[.25 .25 .45 .45])
%subplot(1,3,1)
hold on
U = plot(bootstrap_output.sample_sizes,...
    type_s_errors.uncorrected,...
    'LineStyle','-',...
    'LineWidth',1,...
    'Color',colormap_correlations(1,:));
UA = plot(bootstrap_output.sample_sizes,...
    type_s_errors.average_uncorrected,...
    'LineStyle','-.',...
    'LineWidth',1.5,...
    'Color',colormap_average(1,:));
UB = plot(bootstrap_output.sample_sizes,...
    type_s_errors.best_unthresholded,...
    'LineStyle','-',...
    'LineWidth',1.5,...
    'Color',colormap_average(1,:));

line([bootstrap_output.sample_sizes(required_sample_sizes.uncorrected),...
    bootstrap_output.sample_sizes(required_sample_sizes.uncorrected)],...
    [0,...
    type_s_errors.average_uncorrected(...
    required_sample_sizes.uncorrected)],...
    'LineStyle','--',...
    'LineWidth',1,...
    'Color',[.2 .2 .2]);

line([eps,...
    bootstrap_output.sample_sizes(required_sample_sizes.uncorrected)],...
    [type_s_errors.average_uncorrected(...
    required_sample_sizes.uncorrected),...
    type_s_errors.average_uncorrected(...
    required_sample_sizes.uncorrected)],...
    'LineStyle','--',...
    'LineWidth',1,...
    'Color',[.2 .2 .2]);

line([eps,...
    bootstrap_output.sample_sizes(required_sample_sizes.uncorrected)],...
    [type_s_errors.best_unthresholded(...
    required_sample_sizes.uncorrected),...
    type_s_errors.best_unthresholded(...
    required_sample_sizes.uncorrected)],...
    'LineStyle','--',...
    'LineWidth',1,...
    'Color',[.2 .2 .2]);

hold off
%title('Uncorrected','FontWeight','Normal')
set(gca, 'XScale', xaxis_type,...
    'XTick',sample_size_ticks,...
    'XTickLabel',strsplit(num2str(sample_size_ticks),' '),...
    'YTick',0:20:100,...
    'TickDir','out',...
    'FontName','Arial',...
    'FontSize',figure_font_size,...
    'YGrid','on');
xlim([min(bootstrap_output.sample_sizes),...
    max(bootstrap_output.sample_sizes)])
ylim([0 100])
box off
axis square
set(gcf,'color',[1 1 1])
ylabel('Opposite Sign (%)')
xlabel('Sample size')
legend([U(1),UA(1), UB(1)],...
    {'single','??', 'best'},'location','northwest','Box','off')
pbaspect([1 1 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT INSERT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes('Position',[.52 .46 .2 .2])
plot(bootstrap_output.sample_sizes,...
    type_s_errors.thresholded_uncorrected,...
    'LineStyle','-',...
    'LineWidth',1,...
    'Color',colormap_correlations(1,:));
hold on
plot(bootstrap_output.sample_sizes,...
    type_s_errors.average_thresholded_uncorrected,...
    'LineStyle','-.',...
    'LineWidth',1.5,...
    'Color',colormap_average(1,:));
plot(bootstrap_output.sample_sizes,...
    type_s_errors.best_thresholded_uncorrected,...
    'LineStyle','-',...
    'LineWidth',1.5,...
    'Color',colormap_average(1,:));
set(gca, 'XScale', xaxis_type,...
    'XTick',insert_sample_size_ticks,...
    'XTickLabel',strsplit(num2str(insert_sample_size_ticks),' '),...
    'YTick',[0,3],...
    'TickDir','out',...
    'FontName','Arial',...
    'FontSize',figure_font_size-4,...
    'YGrid','off');
xlim([min(bootstrap_output.sample_sizes),...
    max(bootstrap_output.sample_sizes)])
ylim([0 3])
%axis square
box off
ylabel('Opposite Sign (%)')
xlabel('Sample size')
hold off

if startsWith(lower(save_figures),'y')
    
    image_filename = strcat(output_folder,'/',...
        strrep(dataset_filename,'.csv','_type_s_uncorrected.png'));
    export_fig(image_filename,'-m10','-nocrop', '-transparent','-silent')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FDR CORRECTED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position', figure_size, 'Name','FDR corrected');
axes('Position',[.25 .25 .45 .45])
%subplot(1,3,2)
hold on
F = plot(bootstrap_output.sample_sizes,...
    type_s_errors.fdr,...
    'LineStyle','-',...
    'LineWidth',1,...
    'Color',colormap_correlations(2,:));
FA = plot(bootstrap_output.sample_sizes,...
    type_s_errors.average_fdr,...
    'LineStyle','-.',...
    'LineWidth',1.5,...
    'Color',colormap_average(2,:));
FB = plot(bootstrap_output.sample_sizes,...
    type_s_errors.best_unthresholded,...
    'LineStyle','-',...
    'LineWidth',1.5,...
    'Color',colormap_average(2,:));

line([bootstrap_output.sample_sizes(required_sample_sizes.fdr),...
    bootstrap_output.sample_sizes(required_sample_sizes.fdr)],...
    [0,...
    type_s_errors.average_fdr(...
    required_sample_sizes.fdr)],...
    'LineStyle','--',...
    'LineWidth',1,...
    'Color',[.2 .2 .2]);

line([eps,...
    bootstrap_output.sample_sizes(required_sample_sizes.fdr)],...
    [type_s_errors.average_fdr(...
    required_sample_sizes.fdr),...
    type_s_errors.average_fdr(...
    required_sample_sizes.fdr)],...
    'LineStyle','--',...
    'LineWidth',1,...
    'Color',[.2 .2 .2]);

line([eps,...
    bootstrap_output.sample_sizes(required_sample_sizes.fdr)],...
    [type_s_errors.best_unthresholded(...
    required_sample_sizes.fdr),...
    type_s_errors.best_unthresholded(...
    required_sample_sizes.fdr)],...
    'LineStyle','--',...
    'LineWidth',1,...
    'Color',[.2 .2 .2]);

hold off
%title('FDR Corrected','FontWeight','Normal')
set(gca, 'XScale', xaxis_type,...
    'XTick',sample_size_ticks,...
    'XTickLabel',strsplit(num2str(sample_size_ticks),' '),...
    'YTick',0:20:100,...
    'TickDir','out',...
    'FontName','Arial',...
    'FontSize',figure_font_size,...
    'YGrid','on');
xlim([min(bootstrap_output.sample_sizes),...
    max(bootstrap_output.sample_sizes)])
ylim([0 100])
box off
axis square
set(gcf,'color',[1 1 1])
ylabel('Opposite Sign (%)')
xlabel('Sample size')
legend([F(1),FA(1),FB(1)],...
    {'single','??','best'},'location','northwest','Box','off')
pbaspect([1 1 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT INSERT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes('Position',[.52 .46 .2 .2])
plot(bootstrap_output.sample_sizes,...
    type_s_errors.thresholded_fdr,...
    'LineStyle','-',...
    'LineWidth',1,...
    'Color',colormap_correlations(2,:));
hold on
plot(bootstrap_output.sample_sizes,...
    type_s_errors.average_thresholded_fdr,...
    'LineStyle','-.',...
    'LineWidth',1.5,...
    'Color',colormap_average(2,:));
plot(bootstrap_output.sample_sizes,...
    type_s_errors.best_thresholded_fdr,...
    'LineStyle','-',...
    'LineWidth',1.5,...
    'Color',colormap_average(2,:));
set(gca, 'XScale', xaxis_type,...
    'XTick',insert_sample_size_ticks,...
    'XTickLabel',strsplit(num2str(insert_sample_size_ticks),' '),...
    'YTick',[0,3],...
    'TickDir','out',...
    'FontName','Arial',...
    'FontSize',figure_font_size-4,...
    'YGrid','off');
xlim([min(bootstrap_output.sample_sizes),...
    max(bootstrap_output.sample_sizes)])
ylim([0 3])
%axis square
box off
ylabel('Opposite Sign (%)')
xlabel('Sample size')
hold off

if startsWith(lower(save_figures),'y')
    
    image_filename = strcat(output_folder,'/',...
        strrep(dataset_filename,'.csv','_type_s_fdr_corrected.png'));
    export_fig(image_filename,'-m10','-nocrop', '-transparent','-silent')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% BONFERRONI CORRECTED %%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position', figure_size, 'Name','Bonferroni corrected');
axes('Position',[.25 .25 .45 .45])
%subplot(1,3,3)
hold on
B = plot(bootstrap_output.sample_sizes,...
    type_s_errors.bonferroni,...
    'LineStyle','-',...
    'LineWidth',1,...
    'Color',colormap_correlations(3,:));
BA = plot(bootstrap_output.sample_sizes,...
    type_s_errors.average_bonferroni,...
    'LineStyle','-.',...
    'LineWidth',1.5,...
    'Color',colormap_average(3,:));
BB = plot(bootstrap_output.sample_sizes,...
    type_s_errors.best_unthresholded,...
    'LineStyle','-',...
    'LineWidth',1.5,...
    'Color',colormap_average(3,:));
%title('Bonferroni Corrected','FontWeight','Normal')

line([bootstrap_output.sample_sizes(required_sample_sizes.bonferroni),...
    bootstrap_output.sample_sizes(required_sample_sizes.bonferroni)],...
    [0,...
    type_s_errors.average_bonferroni(...
    required_sample_sizes.bonferroni)],...
    'LineStyle','--',...
    'LineWidth',1,...
    'Color',[.2 .2 .2]);

line([eps,...
    bootstrap_output.sample_sizes(required_sample_sizes.bonferroni)],...
    [type_s_errors.average_bonferroni(...
    required_sample_sizes.bonferroni),...
    type_s_errors.average_bonferroni(...
    required_sample_sizes.bonferroni)],...
    'LineStyle','--',...
    'LineWidth',1,...
    'Color',[.2 .2 .2]);

line([eps,...
    bootstrap_output.sample_sizes(required_sample_sizes.bonferroni)],...
    [type_s_errors.best_unthresholded(...
    required_sample_sizes.bonferroni),...
    type_s_errors.best_unthresholded(...
    required_sample_sizes.bonferroni)],...
    'LineStyle','--',...
    'LineWidth',1,...
    'Color',[.2 .2 .2]);

hold off
set(gca, 'XScale', xaxis_type,...
    'XTick',sample_size_ticks,...
    'XTickLabel',strsplit(num2str(sample_size_ticks),' '),...
    'YTick',0:20:100,...
    'TickDir','out',...
    'FontName','Arial',...
    'FontSize',figure_font_size,...
    'YGrid','on');
xlim([min(bootstrap_output.sample_sizes),...
    max(bootstrap_output.sample_sizes)])
ylim([0 100])
box off
axis square
set(gcf,'color',[1 1 1])
ylabel('Opposite Sign (%)')
xlabel('Sample size')
legend([B(1),BA(1), BB(1)],...
    {'single','??','best'},'location','northwest','Box','off')
pbaspect([1 1 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT INSERT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes('Position',[.52 .46 .2 .2])
plot(bootstrap_output.sample_sizes,...
    type_s_errors.thresholded_bonferroni,...
    'LineStyle','-',...
    'LineWidth',1,...
    'Color',colormap_correlations(3,:));
hold on
plot(bootstrap_output.sample_sizes,...
    type_s_errors.average_thresholded_bonferroni,...
    'LineStyle','-.',...
    'LineWidth',1.5,...
    'Color',colormap_average(3,:));
plot(bootstrap_output.sample_sizes,...
    type_s_errors.best_thresholded_bonferroni,...
    'LineStyle','-',...
    'LineWidth',1.5,...
    'Color',colormap_average(3,:));
set(gca, 'XScale', xaxis_type,...
    'XTick',insert_sample_size_ticks,...
    'XTickLabel',strsplit(num2str(insert_sample_size_ticks),' '),...
    'YTick',[0,3],...
    'TickDir','out',...
    'FontName','Arial',...
    'FontSize',figure_font_size-4,...
    'YGrid','off');
xlim([min(bootstrap_output.sample_sizes),...
    max(bootstrap_output.sample_sizes)])
ylim([0 3])
%axis square
box off
ylabel('Opposite Sign (%)')
xlabel('Sample size')
hold off

% sgtitle(figure_title,'FontWeight','normal',...
%     'FontName','Arial',...
%     'FontSize',figure_font_size+1)

if startsWith(lower(save_figures),'y')
    
    image_filename = strcat(output_folder,'/',...
        strrep(dataset_filename,'.csv','_type_s_bonferroni_corrected.png'));
    export_fig(image_filename,'-m10','-nocrop', '-transparent','-silent')
    
end

end
