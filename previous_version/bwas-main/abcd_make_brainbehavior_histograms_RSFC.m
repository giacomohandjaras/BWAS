%% Make brain behavior histograms (edge, net, component level) for NIH tlbx and cbcl

% make data tables and create sorted matrices and behaviors
abcd_make_big_data_table_FD_08
abcd_make_brainbehavior_variables

% Extract 1st 100 brain components for component analysis below
brain_comp = CorrV(:,1:100);
partitionidx(end + 1) = 333;
% get the edgewise correlations from all factors (ie, factors related to
% cognition)

%% Main
factors = allfactors(:,16);
Edgewise_brbx_correlations_cog = make_edgewise_brbx_RSFC_histogram(allmats_ordered,factors);
Network_brbx_correlations_cog = make_network_brbx_RSFC_histogram(allmats_ordered,factors,partitionidx);
Component_brbx_correlations_cog = make_component_brbx_RSFC_histogram(brain_comp,factors);

% Psychpathology
factors = allfactors(:,30);
Edgewise_brbx_correlations_psy = make_edgewise_brbx_RSFC_histogram(allmats_ordered,factors);
Network_brbx_correlations_psy = make_network_brbx_RSFC_histogram(allmats_ordered,factors,partitionidx);
Component_brbx_correlations_psy = make_component_brbx_RSFC_histogram(brain_comp,factors);

%% Make histogram for Cognition - Main 
fa = .9; % face alpha
ea = 0; % edge alpha
% edgewise
figure;
h1 = histogram(Edgewise_brbx_correlations_cog(:));
h1.Normalization = 'probability';
h1.BinWidth = 0.01;
h1.FaceColor = [199 233 192]./255;
h1.EdgeColor = [199 233 192]./255;
h1.FaceAlpha = fa;
h1.EdgeAlpha = ea;
hold on
% component
h2 = histogram(Component_brbx_correlations_cog(:));
h2.Normalization = 'probability';
h2.BinWidth = 0.01;
h2.FaceColor = [116 196 118]./255;
h2.EdgeColor = [116 196 118]./255;
h2.FaceAlpha = fa;
h2.EdgeAlpha = ea;
% network 
netcorrs = [];
for x = 1:size(Network_brbx_correlations_cog,1)
    for y = x:size(Network_brbx_correlations_cog,2)
        thesecorrs = squeeze(Network_brbx_correlations_cog(x,y,:));
        netcorrs = [netcorrs ; thesecorrs];
    end
end
h3 = histogram(netcorrs);
h3.Normalization = 'probability';
h3.BinWidth = 0.01;
h3.FaceColor = [35 139 69]./255;
h3.EdgeColor = [35 139 69]./255;
h3.FaceAlpha = fa;
h3.EdgeAlpha = ea;

% Legend and axes
legend({'Edge';'Component';'Network'})
%title('RSFC NIH Toolbox Correlations')
%xlabel('Correlation')
%ylabel('Normalized Counts')
xlim([-.2 .2])
ylim([0 .2])
set(gca,'FontSize',16,'FontWeight','bold','FontName','arial')
box('off')
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/RSFCCognitionCorrs_main','tiffn')

%% Make histogram for Clinical vars - Main 
fa = .9; % face alpha
ea = 0; % edge alpha
% edgewise
figure;
h1 = histogram(Edgewise_brbx_correlations_psy(:));
h1.Normalization = 'probability';
h1.BinWidth = 0.01;
h1.FaceColor = [218 218 235]./255;
h1.EdgeColor = [218 218 235]./255;
h1.FaceAlpha = fa;
h1.EdgeAlpha = ea;
hold on
% component
h2 = histogram(Component_brbx_correlations_psy(:));
h2.Normalization = 'probability';
h2.BinWidth = 0.01;
h2.FaceColor = [158 154 200]./255;
h2.EdgeColor = [158 154 200]./255;
h2.FaceAlpha = fa;
h2.EdgeAlpha = ea;
% network 
netcorrs = [];
for x = 1:size(Network_brbx_correlations_psy,1)
    for y = x:size(Network_brbx_correlations_psy,2)
        thesecorrs = squeeze(Network_brbx_correlations_psy(x,y,:));
        netcorrs = [netcorrs ; thesecorrs];
    end
end
h3 = histogram(netcorrs);
h3.Normalization = 'probability';
h3.BinWidth = 0.01;
h3.FaceColor = [106 81 163]./255;
h3.EdgeColor = [106 81 163]./255;
h3.FaceAlpha = fa;
h3.EdgeAlpha = ea;

% Legend and axes
legend({'Edge';'Component';'Network'})
%title('RSFC CBCL Correlations')
%xlabel('Correlation')
%ylabel('Normalized Counts')
xlim([-.2 .2])
set(gca,'FontSize',16,'FontWeight','bold','FontName','arial')
box('off')
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/RSFCClinicalCorrs_main','tiffn')





%% Supplemental 

% Cognition 
factors = allfactors(:,[7:13]);
Edgewise_brbx_correlations_cog = make_edgewise_brbx_RSFC_histogram(allmats_ordered,factors);
Network_brbx_correlations_cog = make_network_brbx_RSFC_histogram(allmats_ordered,factors,partitionidx);
Component_brbx_correlations_cog = make_component_brbx_RSFC_histogram(brain_comp,factors);

% Psychpathology
factors = allfactors(:,20:27);
Edgewise_brbx_correlations_psy = make_edgewise_brbx_RSFC_histogram(allmats_ordered,factors);
Network_brbx_correlations_psy = make_network_brbx_RSFC_histogram(allmats_ordered,factors,partitionidx);
Component_brbx_correlations_psy = make_component_brbx_RSFC_histogram(brain_comp,factors);

%% Make histogram for Cognition - Supplement 
fa = .9; % face alpha
ea = 0; % edge alpha
% edgewise
figure;
h1 = histogram(Edgewise_brbx_correlations_cog(:));
h1.Normalization = 'probability';
h1.BinWidth = 0.01;
h1.FaceColor = [199 233 192]./255;
h1.EdgeColor = [199 233 192]./255;
h1.FaceAlpha = fa;
h1.EdgeAlpha = ea;
hold on
% component
h2 = histogram(Component_brbx_correlations_cog(:));
h2.Normalization = 'probability';
h2.BinWidth = 0.01;
h2.FaceColor = [116 196 118]./255;
h2.EdgeColor = [116 196 118]./255;
h2.FaceAlpha = fa;
h2.EdgeAlpha = ea;
% network 
netcorrs = [];
for x = 1:size(Network_brbx_correlations_cog,1)
    for y = x:size(Network_brbx_correlations_cog,2)
        thesecorrs = squeeze(Network_brbx_correlations_cog(x,y,:));
        netcorrs = [netcorrs ; thesecorrs];
    end
end
h3 = histogram(netcorrs);
h3.Normalization = 'probability';
h3.BinWidth = 0.01;
h3.FaceColor = [35 139 69]./255;
h3.EdgeColor = [35 139 69]./255;
h3.FaceAlpha = fa;
h3.EdgeAlpha = ea;

% Legend and axes
legend({'Edge';'Component';'Network'})
title('RSFC NIH Toolbox Correlations')
xlabel('Correlation')
ylabel('Normalized Counts')
xlim([-.2 .2])
set(gca,'FontSize',16,'FontWeight','bold','FontName','arial')
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/BrainCognitionCorrs_supp','tiffn')

%% Make histogram for Clinical vars - Supplement 
fa = .9; % face alpha
ea = 0; % edge alpha
% edgewise
figure;
h1 = histogram(Edgewise_brbx_correlations_psy(:));
h1.Normalization = 'probability';
h1.BinWidth = 0.01;
h1.FaceColor = [218 218 235]./255;
h1.EdgeColor = [218 218 235]./255;
h1.FaceAlpha = fa;
h1.EdgeAlpha = ea;
hold on
% component
h2 = histogram(Component_brbx_correlations_psy(:));
h2.Normalization = 'probability';
h2.BinWidth = 0.01;
h2.FaceColor = [158 154 200]./255;
h2.EdgeColor = [158 154 200]./255;
h2.FaceAlpha = fa;
h2.EdgeAlpha = ea;
% network 
netcorrs = [];
for x = 1:size(Network_brbx_correlations_psy,1)
    for y = x:size(Network_brbx_correlations_psy,2)
        thesecorrs = squeeze(Network_brbx_correlations_psy(x,y,:));
        netcorrs = [netcorrs ; thesecorrs];
    end
end
h3 = histogram(netcorrs);
h3.Normalization = 'probability';
h3.BinWidth = 0.01;
h3.FaceColor = [106 81 163]./255;
h3.EdgeColor = [106 81 163]./255;
h3.FaceAlpha = fa;
h3.EdgeAlpha = ea;

% Legend and axes
legend({'Edge';'Component';'Network'})
title('RSFC CBCL Correlations')
xlabel('Correlation')
ylabel('Normalized Counts')
xlim([-.2 .2])
set(gca,'FontSize',16,'FontWeight','bold','FontName','arial')
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/BrainCBCLCorrs_supp','tiffn')

