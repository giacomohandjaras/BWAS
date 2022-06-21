%% Cortical Thickness analysis for Figure 1
abcd_make_big_data_table_FD_08

abcd_make_brainbehavior_variables

%% Extract correlation matrices ROI level and vertex level 
% ROI matrices from the discovery and replication sets 
thickness = cell2mat(MainTable_g1.ROIThickness);
thickness_disc = reshape(thickness,333,size(MainTable_g1,1));

thickness = cell2mat(MainTable_g2.ROIThickness);
thickness_rep = reshape(thickness,333,size(MainTable_g2,1));
clear thickness

thickness_disc = thickness_disc(:,MainTable_g1.FrameTotal >= 600 & MainTable_g1.MissingRestSubjects == 0 & MainTable_g1.Badtmaskidx == 0);
% replication set
thickness_rep = thickness_rep(:,MainTable_g2.FrameTotal >= 600 & MainTable_g2.MissingRestSubjects == 0 & MainTable_g2.Badtmaskidx == 0);
allmats_ordered = cat(2,thickness_disc(roi_order(1:333),:),thickness_rep(roi_order(1:333),:));

% Index missing subject 
Missingthickness = zeros(size(allmats_ordered,2),1);
for s = 1:length(Missingthickness)
   thissub = allmats_ordered(:,s);
   if sum(isnan(thissub)) > 0
       Missingthickness(s,1) = 1;
   else
       Missingthickness(s,1) = 0;
   end
end
% Remove any missing subjects from further analysis 
allmats_ordered(:,Missingthickness == 1) = [];
allmats_ordered = allmats_ordered';
allfactors(Missingthickness == 1,:) = [];

% Vertex level 
thickness = cell2mat(MainTable_g1.VertexThickness);
thickness_disc = reshape(thickness,59412,size(MainTable_g1,1));

thickness = cell2mat(MainTable_g2.VertexThickness);
thickness_rep = reshape(thickness,59412,size(MainTable_g2,1));
clear thickness

thickness_disc = thickness_disc(:,MainTable_g1.FrameTotal >= 600 & MainTable_g1.MissingRestSubjects == 0 & MainTable_g1.Badtmaskidx == 0);
% replication set
thickness_rep = thickness_rep(:,MainTable_g2.FrameTotal >= 600 & MainTable_g2.MissingRestSubjects == 0 & MainTable_g2.Badtmaskidx == 0);
allmats_ordered_v = cat(2,thickness_disc,thickness_rep);

% Remove any missing subjects from further analysis 
allmats_ordered_v(:,Missingthickness == 1) = [];
allmats_ordered_v = allmats_ordered_v';

%% Histograms (roi level, network level, component level) Main
% Cog abil 
Vertex_brbx_correlations_cog = make_edge_brbx_structural_histogram(allmats_ordered_v,allfactors(:,16));
ROI_brbx_correlations_cog = make_edge_brbx_structural_histogram(allmats_ordered,allfactors(:,16));
Network_brbx_correlations_cog = make_network_brbx_structural_histogram(allmats_ordered,allfactors(:,16),partitionidx);
% pfactor
Vertex_brbx_correlations_psy = make_edge_brbx_structural_histogram(allmats_ordered_v,allfactors(:,30));
ROI_brbx_correlations_psy = make_edge_brbx_structural_histogram(allmats_ordered,allfactors(:,30));
Network_brbx_correlations_psy = make_network_brbx_structural_histogram(allmats_ordered,allfactors(:,30),partitionidx);

%% Make histogram for Cognition - Main 
fa = .9; % face alpha
ea = 0; % edge alpha

figure;
% Vertex 
h1 = histogram(Vertex_brbx_correlations_cog(:));
h1.Normalization = 'probability';
h1.BinWidth = 0.01;
h1.FaceColor = [199 233 192]./255;
h1.EdgeColor = [199 233 192]./255;
h1.FaceAlpha = 1;
h1.EdgeAlpha = ea;
hold on 
% ROI
h2 = histogram(ROI_brbx_correlations_cog(:));
h2.Normalization = 'probability';
h2.BinWidth = 0.01;
h2.FaceColor = [116 196 118]./255;
h2.EdgeColor = [116 196 118]./255;
h2.FaceAlpha = fa;
h2.EdgeAlpha = ea;
% network 
h3 = histogram(Network_brbx_correlations_cog);
h3.Normalization = 'probability';
h3.BinWidth = 0.01;
h3.FaceColor = [35 139 69]./255;
h3.EdgeColor = [35 139 69]./255;
h3.FaceAlpha = fa;
h3.EdgeAlpha = ea;

% Legend and axes
legend({'Vertex';'ROI';'Network'})
%title('Thickness IQ Correlations')
%xlabel('Correlation')
%ylabel('Normalized Counts')
xlim([-.2 .2])
box('off')
set(gca,'FontSize',16,'FontWeight','bold','FontName','arial')
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/CT_CognitionCorrs_main','tiffn')

%% Make histogram for Clinical vars - Main 
fa = .9; % face alpha
ea = 0; % edge alpha
% edgewise
figure;
h1 = histogram(Vertex_brbx_correlations_psy);
h1.Normalization = 'probability';
h1.BinWidth = 0.01;
h1.FaceColor = [218 218 235]./255;
h1.EdgeColor = [218 218 235]./255;
h1.FaceAlpha = fa;
h1.EdgeAlpha = ea;
hold on
% component
h2 = histogram(ROI_brbx_correlations_psy);
h2.Normalization = 'probability';
h2.BinWidth = 0.01;
h2.FaceColor = [158 154 200]./255;
h2.EdgeColor = [158 154 200]./255;
h2.FaceAlpha = fa;
h2.EdgeAlpha = ea;
% network 
h3 = histogram(Network_brbx_correlations_psy);
h3.Normalization = 'probability';
h3.BinWidth = 0.01;
h3.FaceColor = [106 81 163]./255;
h3.EdgeColor = [106 81 163]./255;
h3.FaceAlpha = fa;
h3.EdgeAlpha = ea;

% Legend and axes
legend({'Vertex';'ROI';'Network'})
%title('Thickness Pfactor Correlations')
%xlabel('Correlation')
%ylabel('Normalized Counts')
xlim([-.2 .2])
box('off')
set(gca,'FontSize',16,'FontWeight','bold','FontName','arial')
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/CT_ClinicalCorrs_main','tiffn')


%% Histograms (roi level, network level, component level)
% IQ 
Vertex_brbx_correlations_cog = make_edge_brbx_structural_histogram(allmats_ordered_v,allfactors(:,[2:11 26]));
ROI_brbx_correlations_cog = make_edge_brbx_structural_histogram(allmats_ordered,allfactors(:,[2:11 26]));
Network_brbx_correlations_cog = make_network_brbx_structural_histogram(allmats_ordered,allfactors(:,[2:11 26]),partitionidx);
% pfactor
Vertex_brbx_correlations_psy = make_edge_brbx_structural_histogram(allmats_ordered_v,allfactors(:,12:21));
ROI_brbx_correlations_psy = make_edge_brbx_structural_histogram(allmats_ordered,allfactors(:,12:21));
Network_brbx_correlations_psy = make_network_brbx_structural_histogram(allmats_ordered,allfactors(:,12:21),partitionidx);
%% Make histogram for Cognition - Supp 
fa = .9; % face alpha
ea = 0; % edge alpha

figure;
% Vertex 
h1 = histogram(Vertex_brbx_correlations_cog(:));
h1.Normalization = 'probability';
h1.BinWidth = 0.01;
h1.FaceColor = [199 233 192]./255;
h1.EdgeColor = [199 233 192]./255;
h1.FaceAlpha = 1;
h1.EdgeAlpha = ea;
hold on 
% ROI
h2 = histogram(ROI_brbx_correlations_cog(:));
h2.Normalization = 'probability';
h2.BinWidth = 0.01;
h2.FaceColor = [116 196 118]./255;
h2.EdgeColor = [116 196 118]./255;
h2.FaceAlpha = fa;
h2.EdgeAlpha = ea;
% network 
h3 = histogram(Network_brbx_correlations_cog);
h3.Normalization = 'probability';
h3.BinWidth = 0.01;
h3.FaceColor = [35 139 69]./255;
h3.EdgeColor = [35 139 69]./255;
h3.FaceAlpha = fa;
h3.EdgeAlpha = ea;

% Legend and axes
legend({'Vertex';'ROI';'Network'})
title('Thickness IQ Correlations')
xlabel('Correlation')
ylabel('Normalized Counts')
xlim([-.2 .2])
set(gca,'FontSize',16,'FontWeight','bold','FontName','arial')
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/BrainStruct_CognitionCorrs_supp','tiffn')

%% Make histogram for Clinical vars - Supp 
fa = .9; % face alpha
ea = 0; % edge alpha
% edgewise
figure;
h1 = histogram(Vertex_brbx_correlations_psy);
h1.Normalization = 'probability';
h1.BinWidth = 0.01;
h1.FaceColor = [218 218 235]./255;
h1.EdgeColor = [218 218 235]./255;
h1.FaceAlpha = fa;
h1.EdgeAlpha = ea;
hold on
% component
h2 = histogram(ROI_brbx_correlations_psy);
h2.Normalization = 'probability';
h2.BinWidth = 0.01;
h2.FaceColor = [158 154 200]./255;
h2.EdgeColor = [158 154 200]./255;
h2.FaceAlpha = fa;
h2.EdgeAlpha = ea;
% network 
h3 = histogram(Network_brbx_correlations_psy);
h3.Normalization = 'probability';
h3.BinWidth = 0.01;
h3.FaceColor = [106 81 163]./255;
h3.EdgeColor = [106 81 163]./255;
h3.FaceAlpha = fa;
h3.EdgeAlpha = ea;

% Legend and axes
legend({'Vertex';'ROI';'Network'})
title('Thickness Pfactor Correlations')
xlabel('Correlation')
ylabel('Normalized Counts')
xlim([-.2 .2])
set(gca,'FontSize',16,'FontWeight','bold','FontName','arial')
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/BrainStruct_ClinicalCorrs_supp','tiffn')

