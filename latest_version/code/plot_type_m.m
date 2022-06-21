function plot_type_m(type_m_errors, bootstrap_output, required_sample_sizes, figure_size, sample_size_ticks, figure_font_size, xaxis_type, colormap_correlations, colormap_average, dataset_filename, save_figures, output_folder)
% coded by Luca Cecchetti & Giacomo Handjaras, IMT-Lucca, Italy
% vers 20220619

%figure_title = strrep(strrep(dataset_filename,'.csv',''),'_',' ');

overall_max_inflation = max(cat(2,...
    max(type_m_errors.inflation_percentage_thresholded_uncorrected(:)),...
    max(type_m_errors.inflation_percentage_thresholded_fdr(:)),...
    max(type_m_errors.inflation_percentage_thresholded_bonferroni(:))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UNCORRECTED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position', figure_size, 'Name','Uncorrected');
axes('Position',[.25 .25 .45 .45])
%subplot(1,3,1)
hold on
U = plot(bootstrap_output.sample_sizes,...
    type_m_errors.inflation_percentage_thresholded_uncorrected,...
    'LineStyle','-',...
    'LineWidth',1,...
    'Color',colormap_correlations(1,:));
UA = plot(bootstrap_output.sample_sizes,...
    type_m_errors.average_inflation_percentage_thresholded_uncorrected,...
    'LineStyle','-.',...
    'LineWidth',1.5,...
    'Color',colormap_average(1,:));
UB = plot(bootstrap_output.sample_sizes,...
    type_m_errors.best_inflation_percentage_thresholded_uncorrected,...
    'LineStyle','-',...
    'LineWidth',1.5,...
    'Color',colormap_average(1,:));
UC = plot(bootstrap_output.sample_sizes,...
    type_m_errors.inflation_percentage_r_crit_uncorrected,...
    'LineStyle',':',...
    'LineWidth',1.5,...
    'Color',[.33 .33 .33]);

line([bootstrap_output.sample_sizes(required_sample_sizes.uncorrected),...
    bootstrap_output.sample_sizes(required_sample_sizes.uncorrected)],...
    [0,...
    type_m_errors.average_inflation_percentage_thresholded_uncorrected(...
    required_sample_sizes.uncorrected)],...
    'LineStyle','--',...
    'LineWidth',1,...
    'Color',[.2 .2 .2]);

line([eps,...
    bootstrap_output.sample_sizes(required_sample_sizes.uncorrected)],...
    [type_m_errors.average_inflation_percentage_thresholded_uncorrected(...
    required_sample_sizes.uncorrected),...
    type_m_errors.average_inflation_percentage_thresholded_uncorrected(...
    required_sample_sizes.uncorrected)],...
    'LineStyle','--',...
    'LineWidth',1,...
    'Color',[.2 .2 .2]);

line([eps,...
    bootstrap_output.sample_sizes(required_sample_sizes.uncorrected)],...
    [type_m_errors.best_inflation_percentage_thresholded_uncorrected(...
    required_sample_sizes.uncorrected),...
    type_m_errors.best_inflation_percentage_thresholded_uncorrected(...
    required_sample_sizes.uncorrected)],...
    'LineStyle','--',...
    'LineWidth',1,...
    'Color',[.2 .2 .2]);
hold off

%title('Uncorrected','FontWeight','Normal')
set(gca, 'XScale', xaxis_type,...
    'XTick',sample_size_ticks,...
    'XTickLabel',strsplit(num2str(sample_size_ticks),' '),...
    'YTick',0:250:overall_max_inflation,...
    'TickDir','out',...
    'FontName','Arial',...
    'FontSize',figure_font_size,...
    'YGrid','on');
xlim([min(bootstrap_output.sample_sizes),...
    max(bootstrap_output.sample_sizes)])
ylim([0 overall_max_inflation])
box off
%axis square
set(gcf,'color',[1 1 1])
ylabel('Inflation (%)')
xlabel('Sample size')
legend([U(1),UA(1),UB(1),UC(1)],...
    {'single','µ','best','crit'},'location','northeast','Box','off')
%pbaspect([1 1 1])

if startsWith(lower(save_figures),'y')
    
    image_filename = strcat(output_folder,'/',...
        strrep(dataset_filename,'.csv','_type_m_uncorrected.png'));
    export_fig(image_filename,'-m10','-nocrop', '-transparent','-silent')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FDR CORRECTED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position', figure_size, 'Name','FDR corrected');
axes('Position',[.25 .25 .45 .45])
%subplot(1,3,2)
hold on
F = plot(bootstrap_output.sample_sizes,...
    type_m_errors.inflation_percentage_thresholded_fdr,...
    'LineStyle','-',...
    'LineWidth',1,...
    'Color',colormap_correlations(2,:));
FA = plot(bootstrap_output.sample_sizes,...
    type_m_errors.average_inflation_percentage_thresholded_fdr,...
    'LineStyle','-.',...
    'LineWidth',1.5,...
    'Color',colormap_average(2,:));
FB = plot(bootstrap_output.sample_sizes,...
    type_m_errors.best_inflation_percentage_thresholded_fdr,...
    'LineStyle','-',...
    'LineWidth',1.5,...
    'Color',colormap_average(2,:));

line([bootstrap_output.sample_sizes(required_sample_sizes.fdr),...
    bootstrap_output.sample_sizes(required_sample_sizes.fdr)],...
    [0,...
    type_m_errors.average_inflation_percentage_thresholded_fdr(...
    required_sample_sizes.fdr)],...
    'LineStyle','--',...
    'LineWidth',1,...
    'Color',[.2 .2 .2]);

line([eps,...
    bootstrap_output.sample_sizes(required_sample_sizes.fdr)],...
    [type_m_errors.average_inflation_percentage_thresholded_fdr(...
    required_sample_sizes.fdr),...
    type_m_errors.average_inflation_percentage_thresholded_fdr(...
    required_sample_sizes.fdr)],...
    'LineStyle','--',...
    'LineWidth',1,...
    'Color',[.2 .2 .2]);

line([eps,...
    bootstrap_output.sample_sizes(required_sample_sizes.fdr)],...
    [type_m_errors.best_inflation_percentage_thresholded_fdr(...
    required_sample_sizes.fdr),...
    type_m_errors.best_inflation_percentage_thresholded_fdr(...
    required_sample_sizes.fdr)],...
    'LineStyle','--',...
    'LineWidth',1,...
    'Color',[.2 .2 .2]);

%title('FDR Corrected','FontWeight','Normal')
set(gca, 'XScale', xaxis_type,...
    'XTick',sample_size_ticks,...
    'XTickLabel',strsplit(num2str(sample_size_ticks),' '),...
    'YTick',0:250:overall_max_inflation,...
    'TickDir','out',...
    'FontName','Arial',...
    'FontSize',figure_font_size,...
    'YGrid','on');
xlim([min(bootstrap_output.sample_sizes),...
    max(bootstrap_output.sample_sizes)])
ylim([0 overall_max_inflation])
box off
%axis square
set(gcf,'color',[1 1 1])
ylabel('Inflation (%)')
xlabel('Sample size')
legend([F(1),FA(1),FB(1)],...
    {'single','µ','best'},'location','northeast','Box','off')
%pbaspect([1 1 1])

if startsWith(lower(save_figures),'y')
    
    image_filename = strcat(output_folder,'/',...
        strrep(dataset_filename,'.csv','_type_m_fdr_corrected.png'));
    export_fig(image_filename,'-m10','-nocrop', '-transparent','-silent')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% BONFERRONI CORRECTED %%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position', figure_size, 'Name','Bonferroni corrected');
axes('Position',[.25 .25 .45 .45])
%subplot(1,3,3)
hold on
B = plot(bootstrap_output.sample_sizes,...
    type_m_errors.inflation_percentage_thresholded_bonferroni,...
    'LineStyle','-',...
    'LineWidth',1,...
    'Color',colormap_correlations(3,:));
BA = plot(bootstrap_output.sample_sizes,...
    type_m_errors.average_inflation_percentage_thresholded_bonferroni,...
    'LineStyle','-.',...
    'LineWidth',1.5,...
    'Color',colormap_average(3,:));
BB = plot(bootstrap_output.sample_sizes,...
    type_m_errors.best_inflation_percentage_thresholded_bonferroni,...
    'LineStyle','-',...
    'LineWidth',1.5,...
    'Color',colormap_average(3,:));
BC = plot(bootstrap_output.sample_sizes,...
    type_m_errors.inflation_percentage_r_crit_bonferroni,...
    'LineStyle',':',...
    'LineWidth',1.5,...
    'Color',[.33 .33 .33]);

line([bootstrap_output.sample_sizes(required_sample_sizes.bonferroni),...
    bootstrap_output.sample_sizes(required_sample_sizes.bonferroni)],...
    [0,...
    type_m_errors.average_inflation_percentage_thresholded_bonferroni(...
    required_sample_sizes.bonferroni)],...
    'LineStyle','--',...
    'LineWidth',1,...
    'Color',[.2 .2 .2]);

line([eps,...
    bootstrap_output.sample_sizes(required_sample_sizes.bonferroni)],...
    [type_m_errors.average_inflation_percentage_thresholded_bonferroni(...
    required_sample_sizes.bonferroni),...
    type_m_errors.average_inflation_percentage_thresholded_bonferroni(...
    required_sample_sizes.bonferroni)],...
    'LineStyle','--',...
    'LineWidth',1,...
    'Color',[.2 .2 .2]);

line([eps,...
    bootstrap_output.sample_sizes(required_sample_sizes.bonferroni)],...
    [type_m_errors.best_inflation_percentage_thresholded_bonferroni(...
    required_sample_sizes.bonferroni),...
    type_m_errors.best_inflation_percentage_thresholded_bonferroni(...
    required_sample_sizes.bonferroni)],...
    'LineStyle','--',...
    'LineWidth',1,...
    'Color',[.2 .2 .2]);

%title('Bonferroni Corrected','FontWeight','Normal')
set(gca, 'XScale', xaxis_type,...
    'XTick',sample_size_ticks,...
    'XTickLabel',strsplit(num2str(sample_size_ticks),' '),...
    'YTick',0:250:overall_max_inflation,...
    'TickDir','out',...
    'FontName','Arial',...
    'FontSize',figure_font_size,...
    'YGrid','on');
xlim([min(bootstrap_output.sample_sizes),...
    max(bootstrap_output.sample_sizes)])
ylim([0 overall_max_inflation])
box off
%axis square
set(gcf,'color',[1 1 1])
ylabel('Inflation (%)')
xlabel('Sample size')
legend([B(1),BA(1),BB(1),BC(1)],...
    {'single','µ','best','crit'},'location','northeast','Box','off')
%pbaspect([1 1 1])

% sgtitle(figure_title,'FontWeight','normal',...
%     'FontName','Arial',...
%     'FontSize',figure_font_size+1)

if startsWith(lower(save_figures),'y')
    
    image_filename = strcat(output_folder,'/',...
        strrep(dataset_filename,'.csv','_type_m_bonferroni_corrected.png'));
    export_fig(image_filename,'-m10','-nocrop', '-transparent','-silent')
    
end

end


