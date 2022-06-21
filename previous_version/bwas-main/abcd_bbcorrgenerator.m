function [component_output,vertex_output,numbersubjects] = abcd_bbcorrgenerator(table)
% input is a text file generated from abcde capture tool
% format is subject ID (column 1) followed by behavioral vars of interest 
% Function reads in pre-existing brain data table and uses innerjoin to
% join with behavioral table 
addpath(genpath('/data/cn/data1/scripts/CIFTI_RELATED/'))
addpath('/data/nil-bluearc/GMT/Scott/ABCD/Scripts/')
addpath(genpath('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/ABCD_pitt_washu/Scripts/'))

% Read subject list - all the rest subjects we have 
RestTable = readtable('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/SubjectLists/total_subjectlist_0708.txt');
RestTable.Properties.VariableNames{1} = 'Subject';
subjects = RestTable.Subject;


%% Append correlation matrices & motion information  
%load('/data/nil-bluearc/GMT/Scott/ABCD/Parcelcorrmats/SubjectParcelCorrMats.mat')
load('/data/nil-bluearc/GMT/Scott/ABCD/Parcelcorrmats/SubjectParcelCorrMats_withsubcort_0.08FD.mat')
%load('/data/nil-bluearc/GMT/Scott/ABCD/Parcelcorrmats/SubjectParcelCorrMats_withsubcort_noFD_0.08FD.mat');
load('/data/nil-bluearc/GMT/Scott/ABCD/FD/SubjectFD_0.08.mat');
load('/data/nil-bluearc/GMT/Scott/ABCD/FD/SubjectFD_Incframes_0.08.mat')
load('/data/nil-bluearc/GMT/Scott/ABCD/FD/SubjectFrameTotal_0.08.mat')
load('/data/nil-bluearc/GMT/Scott/ABCD/FD/Badtmaskidx.mat')
load('/data/nil-bluearc/GMT/Scott/ABCD/Parcelcorrmats/MissingSubIdx.mat')
load('/data/nil-bluearc/GMT/Scott/ABCD/FD/FDtimecourse.mat');
load('/data/nil-bluearc/GMT/Scott/ABCD/Parcelcorrmats/SplitHalfReliability_rest.mat')
load('/data/nil-bluearc/GMT/Scott/ABCD/CorticalThickness/CorticalThickness.mat');
load('/data/nil-bluearc/GMT/Scott/ABCD/CorticalThickness/MissingSubIdx_CortThickness.mat');  
load('/data/nil-bluearc/GMT/Scott/ABCD/CorticalThickness/ParcelThickness.mat');
IncompleteRest = zeros(size(subjects,1),1);

for i = 1:size(RestTable,1)
    Thismat = Corrmat(:,:,i);
    %nofdmat = Corrmat_noFD(:,:,i);
    VertexThickness{i,1} = CorticalThickness(:,i);
    ROIThickness{i,1} = ParcelThickness(:,i);
    % check for nans 
    if sum(any(isnan(Thismat))) > 0
        IncompleteRest(i,1) = 1;
    else
        Corrmatcell{i,1} = Thismat;
        %Corrmatcell_noFD{i,1} = nofdmat;
    end

end

RestTable(:,2) = Corrmatcell;
%RestTable(:,3) = Corrmatcell_noFD;
RestTable(:,3) = num2cell(SubjectFD);
RestTable(:,4) = num2cell(SubjectFD_Incframes);
RestTable(:,5) = num2cell(SubjectFrameTotal);
RestTable(:,6) = num2cell(MissingSubIdx);
RestTable(:,7) = num2cell(Badtmaskidx);
RestTable(:,8) = num2cell(FDtimecourse);
RestTable(:,9) = num2cell(IncompleteRest);
RestTable(:,10) = num2cell(Reliability);
RestTable(:,11) = VertexThickness;
RestTable(:,12) = ROIThickness;
RestTable(:,13) = num2cell(MissingSubIdx_CortThick);
RestTable.Properties.VariableNames{2} = 'Corrmats';
%RestTable.Properties.VariableNames{3} = 'Corrmats_noFD';
RestTable.Properties.VariableNames{3} = 'FD';
RestTable.Properties.VariableNames{4} = 'FD_Incframes';
RestTable.Properties.VariableNames{5} = 'FrameTotal';
RestTable.Properties.VariableNames{6} = 'MissingRestSubjects';
RestTable.Properties.VariableNames{7} = 'Badtmaskidx';
RestTable.Properties.VariableNames{8} = 'FDtimecourse';
RestTable.Properties.VariableNames{9} = 'IncompleteRest';
RestTable.Properties.VariableNames{10} = 'Reliability';
RestTable.Properties.VariableNames{11} = 'VertexThickness';
RestTable.Properties.VariableNames{12} = 'ROIThickness';
RestTable.Properties.VariableNames{13} = 'MissingThickness';

%% Remove subjects with incomplete rest data 
RestTable(RestTable.IncompleteRest == 1 | RestTable.FrameTotal < 600 | RestTable.MissingRestSubjects == 1 | RestTable.MissingThickness == 1 | RestTable.Badtmaskidx == 1,:) = [];


%% Read table from abcde output 
behavioraltable = readtable(tablefile);
behavioraltable.x_subjectkey=arrayfun(@(s)(extractAfter(s(1),'NDAR_')),behavioraltable.x_subjectkey);
behavioraltable.Properties.VariableNames{1} = 'Subject';
%cu = load('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/CUtable.mat'); 
%cu = cu.CUtable;
%behavioraltable = innerjoin(behavioraltable,cu,'Keys','Subject');
%behavioraltable = movevars(behavioraltable,'interview_age','Before','height');
behavioraltable.income(behavioraltable.income==999 | behavioraltable.income==777) = nan;
behavioraltable.numpeeps(behavioraltable.numpeeps>10) = 11;
income2needs = behavioraltable.income./behavioraltable.numpeeps;
income2needs(isinf(income2needs)) = nan;
behavioraltable = addvars(behavioraltable,income2needs,'After','income');
behavioraltable.numpeeps = [];
behavioraltable.income = [];

% Remove missing data 
id = {'NaN' '' 999 777 NaN Inf '999' '777'};
TF = ismissing(behavioraltable,id);

TFs = sum(TF')';

% Remove subjects with missing data 
behavioraltable(TFs>0,:) = [];

% Join with resting state data 
RestTable = innerjoin(RestTable,behavioraltable,'Keys','Subject');

clear Corrmat nofdmat ParcelThickness CorticalThickness 
vars = table2array(behavioraltable(:,2:end))';
% make correlation matrix of bahavioral variables (sans rest)
withinbehaviorcorrs = zeros(size(vars,1));
for i = 1:size(vars,1)
    for j = i+1:size(vars,1)
            withinbehaviorcorrs(i,j) = corr(vars(i,:)',vars(j,:)','Type','Pearson');
    end
end
withinbehaviorcorrs = withinbehaviorcorrs + withinbehaviorcorrs';
cmap = [];
cmap(3,:) = [177 0 0]./255;
cmap(2,:) = [255 255 255]./255;
cmap(1,:) = [12 44 132]./255;
[x,y]=meshgrid([1:size(cmap,1)],1:200);
cmap = interp2(x([1,100,200],:),y([1,100,200],:),cmap,x,y);
figure; imagesc(withinbehaviorcorrs);colormap(cmap);colorbar;caxis([-1 1])
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/All_bahevioral_corrs_supp','tiffn')

%% Convert 3D matrices isnto 2D (edge x subject) and also extract components vis svd 
rmats = cell2mat(RestTable.Corrmats);
rmats = reshape(rmats,size(rmats,2),size(RestTable,1),size(rmats,2));
rmats = permute(rmats,[1,3,2]);

% Reorder matrices according to Gordon/Laumann network partitions
reorder_gordon_laumann_parcels
% Get ROI order 
roi_order = NetworksOrdered(:,1);
%Reorder
roi_order(334:size(rmats,1)) = 334:394;

rmats_ordered = rmats(roi_order,roi_order,:);

% convert 3d matrices into 2d 
uidx = find(triu(rmats(:,:,1),1));
allmats_2d = zeros(length(uidx),size(rmats,3));
for s = 1:size(allmats_2d,2)
    thisr = rmats(:,:,s);
    allmats_2d(:,s) = thisr(uidx);
end

% Behavioral variables that have rest data and cort thickness data  
vars = table2array(RestTable(:,14:end));
allmats_2d = allmats_2d';

% Initialize output 
component_output = cell(size(vars,2),1);
vertex_output = cell(size(vars,2),1);

[~,braincomponents,~,~,explained] = pca(allmats_2d);

z = parpool(15);
parfor v = 1:size(vars,2)
    disp(['On variable ' num2str(v)])
    
    %component wise corrs 
    thesecorrs = zeros(size(braincomponents,1),1);
    for c = 1:size(braincomponents,2)
         thesecorrs(c,1) = corr(braincomponents(:,c),vars(:,v));    
    end
    component_output{v} = thesecorrs;
    
    % vertex wise corrs
    thesecorrs = zeros(size(allmats_2d,2),1);
    for e = 1:size(allmats_2d,2)
         thesecorrs(e,1) = corr(allmats_2d(:,e),vars(:,v));  
    end
    vertex_output{v} = thesecorrs;
    
end
delete(z)

% Networks 
partitionidx(end+1) = 333;
network_output = make_network_brbx_RSFC_histogram(rmats_ordered,vars,partitionidx);

% Resize to 2d 
corrs2d = [];
uidx = find(triu(network_output(:,:,1)));
for i = 1:size(vars,2)
    thesecorrs = network_output(:,:,i);
    corrs2d(:,i) = thesecorrs(uidx);
end
network_output = corrs2d;

% Save into structure 
AllBBcorrsRSFC.component_output = component_output;
AllBBcorrsRSFC.vertex_output = vertex_output;
AllBBcorrsRSFC.network_output = network_output;
save('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/AllBBcorrsRSFC.mat','AllBBcorrsRSFC')


%% Make all subplots 
load('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/AllBBcorrsRSFC.mat','AllBBcorrsRSFC')

% Make X x 3 color palette
colors = [];
% Demographics 
cmap = [];
cmap(3,:) = [8 69 148]./255;
cmap(2,:) = [66 146 198]./255;
cmap(1,:) = [158 202 225]./255;
[x,y]=meshgrid([1:size(cmap,1)],1:6);
cmap = interp2(x([1,3,6],:),y([1,3,6],:),cmap,x,y);

colors = cmap;

% Cognition 
cmap = [];
cmap(3,:) = [0 90 50]./255;
cmap(2,:) = [65 171 93]./255;
cmap(1,:) = [161 217 155]./255;
[x,y]=meshgrid([1:size(cmap,1)],1:13);
cmap = interp2(x([1,7,13],:),y([1,7,13],:),cmap,x,y);

colors = [colors;cmap];

% Mental health 
cmap = [];
cmap(3,:) = [74 20 134]./255;
cmap(2,:) = [128 125 186]./255;
cmap(1,:) = [188 189 200]./255;
[x,y]=meshgrid([1:size(cmap,1)],1:22);
cmap = interp2(x([1,11,22],:),y([1,11,22],:),cmap,x,y);

colors = [colors;cmap];

% Plot 
fa = .9;
ea = 0;
figure; 
for v = 1:size(AllBBcorrsRSFC.component_output,1)
    corrs = [];
    corrs = [cell2mat(AllBBcorrsRSFC.component_output(v)) ; cell2mat(AllBBcorrsRSFC.vertex_output(v))]; 
    subplot(7,7,v);
    h1 = histogram(corrs);
    h1.Normalization = 'probability';
    h1.BinWidth = 0.01;
    h1.FaceColor = colors(v,:);
    h1.EdgeColor = colors(v,:);
    h1.FaceAlpha = fa;
    h1.EdgeAlpha = ea;
    xlim([-.15 .15])
    set(gca,'XTick',[],'YTick',[])
    set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
    box('off')
end
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/RSFC_allbbcorrs_supp','tiffn')

%% structural data 
% Extract correlation matrices ROI level and vertex level 

% Vertex level 
thickness = cell2mat(RestTable.VertexThickness);
allmats_ordered_v  = reshape(thickness,59412,size(RestTable,1));
clear thickness

Missingthickness = zeros(size(allmats_ordered_v,2),1);
for s = 1:length(Missingthickness)
   thissub = allmats_ordered_v(:,s);
   if sum(isnan(thissub)) > 0
       Missingthickness(s,1) = 1;
   else
       Missingthickness(s,1) = 0;
   end
end

% Remove any missing subjects from further analysis 
allmats_ordered_v(:,Missingthickness == 1) = [];
allmats_ordered_v = allmats_ordered_v';
% Behavioral variables that have rest data and cort thickness data  
vars = table2array(RestTable(:,14:end));
vars(Missingthickness == 1,:) = [];

% Initialize output 
component_output = cell(size(vars,2),1);
vertex_output = cell(size(vars,2),1);

[~,braincomponents,~,~,explained] = pca(allmats_ordered_v);

z = parpool(15);
parfor v = 1:size(vars,2)
    
    tic;
    disp(['On variable ' num2str(v)])
    % pca / component wise corrs 
    thesecorrs = zeros(size(braincomponents,2),1);
    for c = 1:size(braincomponents,2)
         thesecorrs(c,1) = corr(braincomponents(:,c),vars(:,v));    
    end
    component_output{v} = thesecorrs;
    
    % vertex wise corrs
    thesecorrs = zeros(size(allmats_ordered_v,2),1);
    for e = 1:size(allmats_ordered_v,2)
         thesecorrs(e,1) = corr(allmats_ordered_v(:,e),vars(:,v));  
    end
    vertex_output{v} = thesecorrs;
    
    toc;   
end
delete(z)

% Save into structure 
AllBBcorrsthickness.component_output = component_output;
AllBBcorrsthickness.vertex_output = vertex_output;
save('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/AllBBcorrsthickness.mat','AllBBcorrsthickness')

%% 
load('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/AllBBcorrsthickness.mat')

% Make X x 3 color palette
colors = [];
% Demographics 
cmap = [];
cmap(3,:) = [8 69 148]./255;
cmap(2,:) = [66 146 198]./255;
cmap(1,:) = [158 202 225]./255;
[x,y]=meshgrid([1:size(cmap,1)],1:6);
cmap = interp2(x([1,3,6],:),y([1,3,6],:),cmap,x,y);

colors = cmap;

% Cognition 
cmap = [];
cmap(3,:) = [0 90 50]./255;
cmap(2,:) = [65 171 93]./255;
cmap(1,:) = [161 217 155]./255;
[x,y]=meshgrid([1:size(cmap,1)],1:13);
cmap = interp2(x([1,7,13],:),y([1,7,13],:),cmap,x,y);

colors = [colors;cmap];

% Mental health 
cmap = [];
cmap(3,:) = [74 20 134]./255;
cmap(2,:) = [128 125 186]./255;
cmap(1,:) = [188 189 200]./255;
[x,y]=meshgrid([1:size(cmap,1)],1:22);
cmap = interp2(x([1,11,22],:),y([1,11,22],:),cmap,x,y);

colors = [colors;cmap];

% Plot 
fa = .9;
ea = 0;
figure; 
for v = 1:size(AllBBcorrsthickness.component_output,1)
    corrs = [];
    corrs = [cell2mat(AllBBcorrsthickness.component_output(v)) ; cell2mat(AllBBcorrsthickness.vertex_output(v))]; 
    subplot(7,7,v);
    h1 = histogram(corrs);
    h1.Normalization = 'probability';
    h1.BinWidth = 0.01;
    h1.FaceColor = colors(v,:);
    h1.EdgeColor = colors(v,:);
    h1.FaceAlpha = fa;
    h1.EdgeAlpha = ea;
    xlim([-.15 .15])
    set(gca,'XTick',[],'YTick',[])
    set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
    box('off')
end
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/CT_allbbcorrs_supp','tiffn')


%% Get 1st percentile of vertex and edge corrs across everything 
load('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/AllBBcorrsRSFC.mat')
load('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/AllBBcorrsthickness.mat')
allcorrs = [];
for v = 1:size(AllBBcorrsRSFC.component_output,1)
    allcorrs = [allcorrs;cell2mat(AllBBcorrsRSFC.vertex_output(v)); cell2mat(AllBBcorrsRSFC.component_output(v)) ; AllBBcorrsRSFC.network_output(:,v)]; 
    %allcorrs = [allcorrs;cell2mat(AllBBcorrsthickness.vertex_output(v)) ; cell2mat(AllBBcorrsthickness.component_output(v))];
end
prctile(abs(allcorrs),99)
min(allcorrs)
max(allcorrs)
nanmedian(abs(allcorrs))


%% Get 1st percentile of IQ and Pfactor across modalities 
load('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/AllBBcorrsRSFC.mat')
load('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/AllBBcorrsthickness.mat')

% Iq = 17 ; across RSFC and thickness
alliqcorrs = [cell2mat(AllBBcorrsthickness.vertex_output(17));cell2mat(AllBBcorrsRSFC.vertex_output(17))]; 
prctile(abs(alliqcorrs),99)
% P = 30  ; across RSFC and thickness
allpfcorrs = [cell2mat(AllBBcorrsthickness.vertex_output(30));cell2mat(AllBBcorrsRSFC.vertex_output(30))];
prctile(abs(allpfcorrs),99)

%RSFC ; across IQ and pfactor 
allrsfccorrs = [cell2mat(AllBBcorrsRSFC.vertex_output(17));cell2mat(AllBBcorrsRSFC.vertex_output(30))]; 
prctile(abs(allrsfccorrs),99)

%Thickness ; across IQ and pfactor 
allrsfccorrs = [cell2mat(AllBBcorrsthickness.vertex_output(17));cell2mat(AllBBcorrsthickness.vertex_output(30))]; 
prctile(abs(allrsfccorrs),99)

end
