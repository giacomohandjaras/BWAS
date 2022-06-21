function plot_correlation_histogram(simulated_data,n_bins,dataset_filename,figure_size,figure_font_size,plot_colormap,save_figures,output_folder)
% coded by Luca Cecchetti & Giacomo Handjaras, IMT-Lucca, Italy
% vers 20220619

%figure_title = strrep(strrep(dataset_filename,'.csv',''),'_',' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ORIGINAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Position', figure_size, 'Name','Actual data');
axes('Position',[.25 .25 .45 .45])
%subplot(1,2,1)
H1 = histogram(simulated_data.original_corr_values,n_bins,...
    'Normalization','probability',...
    'FaceColor',plot_colormap(1,:),...
    'FaceAlpha',0.5,...
    'EdgeColor','none');
xlim([-max(abs(simulated_data.original_corr_values))-H1.BinWidth,...
    max(abs(simulated_data.original_corr_values))+H1.BinWidth]);
ylim([0 max(H1.BinCounts./sum(H1.BinCounts))])
xlabel('Correlation (r)')
ylabel('Occurrence (p)')
set(gca,'FontName','Arial',...
    'FontSize',figure_font_size,...
    'YTick',0:0.1:max(H1.BinCounts./sum(H1.BinCounts)),...
    'TickDir','out')
box off
axis square
set(gcf,'Color',[1 1 1])
%title('Actual data','FontWeight','normal')
pbaspect([1 1 1])

if startsWith(lower(save_figures),'y')
    
    image_filename = strcat(output_folder,'/',...
        strrep(dataset_filename,'.csv','_actual_distribution.png'));
    export_fig(image_filename,'-m10','-nocrop', '-transparent','-silent')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position', figure_size, 'Name','Simulated data');
axes('Position',[.25 .25 .45 .45])
%subplot(1,2,2)
H2 = histogram(simulated_data.simulated_corr_values,n_bins,...
    'Normalization','probability',...
    'FaceColor',plot_colormap(1,:),...
    'FaceAlpha',0.5,...
    'EdgeColor','none');
xlim([-max(abs(simulated_data.original_corr_values))-H1.BinWidth,...
    max(abs(simulated_data.original_corr_values))+H1.BinWidth]);
ylim([0 max(H1.BinCounts./sum(H1.BinCounts))])
xlabel('Correlation (r)')
ylabel('Occurrence (p)')
set(gca,'FontName','Arial',...
    'FontSize',figure_font_size,...
    'YTick',0:0.1:max(H1.BinCounts./sum(H1.BinCounts)),...
    'TickDir','out')
box off
axis square
set(gcf,'Color',[1 1 1])
%title('Simulated data','FontWeight','normal')
pbaspect([1 1 1])

% sgtitle(figure_title,'FontWeight','normal',...
%     'FontName','Arial',...
%     'FontSize',figure_font_size)

if startsWith(lower(save_figures),'y')
    
    image_filename = strcat(output_folder,'/',...
        strrep(dataset_filename,'.csv','_simulated_distribution.png'));
    export_fig(image_filename,'-m10','-nocrop', '-transparent','-silent')
    
end


end
