%% Make vibration plot for FIgure 1B, !C and 1E, 1F 
% Includes both RSFC and thickness 


%% Figure 1f sampling variability RSFC/IQ
addpath('/data/nil-bluearc/GMT/Scott/MSC_Subcortical/Scripts')
abcd_make_big_data_table_FD_08

abcd_make_brainbehavior_variables

%% Cognition 
% Run all brain behavior correlations
factors = allfactors(:,16); % run across these factors 
FullSampleCorrs = zeros(size(allmats_2d,2),size(factors,2));
for x = 1:size(allmats_2d,2)
    for f = 1:size(factors,2)
        FullSampleCorrs(x,f) = corr(squeeze(allmats_2d(:,x)),factors(:,f));
    end
end

% Get max abs correlation
[a,idx] = max(abs(FullSampleCorrs));

% % Get top 1% ABSOLUTE value corrs
% mat = FullSampleCorrs(:,:,1);
% Input = zeros(size(allmats_ordered,2),floor(size(mat,1)/100),size(FullSampleCorrs,2)); % Initiate 2D matrices
% for v = 1:size(FullSampleCorrs,2)
%     Thisvar = FullSampleCorrs(:,:,v);
%     [a b] = sort(abs(Thisvar),'descend');
%     % Get top 1% vertices across all subs
%     Input(:,:,v) = allmats_ordered(b(1:floor(size(mat,1)/100)),:)';
% end
% rmats = Input; clear Input

% Set input params  
binsize = ([25 33 50 70 100 135 200 265 375 525 725 1000 1430 2000 2800 3928]);

iter = 1000; % how many times to draw with replacement at each interval 

% Get edge correlations resampling incrementally larger samples 
f = 16;
[Corrs,~] = abcd_edgewise_correlation_iterative_reliability_single_factor(allmats_2d(:,idx),allfactors(:,f),binsize,iter);
Corrs = squeeze(Corrs);
prctileplotting(Corrs,[116 196 118],'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/SV_rsfcIQ_CIs')


%% Clinical 
factors = allfactors(:,30); % run across these factors 
FullSampleCorrs = zeros(size(allmats_2d,2),size(factors,2));
for x = 1:size(allmats_2d,2)
    for f = 1:size(factors,2)
        FullSampleCorrs(x,f) = corr(squeeze(allmats_2d(:,x)),factors(:,f));
    end
end

% Get max abs correlation
[a,idx] = max(abs(FullSampleCorrs));

% % Get top 1% ABSOLUTE value corrs
% mat = FullSampleCorrs(:,:,1);
% Input = zeros(size(allmats_ordered,2),floor(size(mat,1)/100),size(FullSampleCorrs,2)); % Initiate 2D matrices
% for v = 1:size(FullSampleCorrs,2)
%     Thisvar = FullSampleCorrs(:,:,v);
%     [a b] = sort(abs(Thisvar),'descend');
%     % Get top 1% vertices across all subs
%     Input(:,:,v) = allmats_ordered(b(1:floor(size(mat,1)/100)),:)';
% end
% rmats = Input; clear Input

% Set input params  
binsize = ([25 33 50 70 100 135 200 265 375 525 725 1000 1430 2000 2800 3928]);
iter = 1000; % how many times to draw with replacement at each interval 

% Get edge correlations resampling incrementally larger samples 
f = 30;
[Corrs,~] = abcd_edgewise_correlation_iterative_reliability_single_factor(allmats_2d(:,idx),allfactors(:,f),binsize,iter);
Corrs = squeeze(Corrs);
prctileplotting(Corrs,[158 154 200],'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/SV_rsfcPF_CIs')
%% Scatter plot (Fig 1C)

abcd_make_big_data_table_FD_08

abcd_make_brainbehavior_variables

% Grab DMN corrs for each subject 
DMNarray = [];
for s = 1:size(allmats_ordered,3)
   if s == 1
       thissub = allmats_ordered(1:41,1:41,s);
       uidx = find(triu(thissub,1));
       DMNarray(:,s) = thissub(uidx);
   else
       thissub = allmats_ordered(1:41,1:41,s);
       DMNarray(:,s) = thissub(uidx);
   end    
end
% Full sample correlation and p value 
mDMN = mean(DMNarray,1)';
[r p] = corr(mDMN,allfactors(:,16));
disp(['Full sample DMN RSFC IQ correlation = ' num2str(r) ' : p value = ' num2str(p)]);
[r p] = corr(mDMN,allfactors(:,30));
disp(['Full sample DMN RSFC p factor correlation = ' num2str(r) ' : p value = ' num2str(p)]);

% Cognition 
% Loop thru 1000 iterations of this correlation, grabbing 25 subs each time
% Save out subject indices at each loop 
for i = 1:1000
    idx = datasample(1:size(allmats_ordered,3),25)';
    iteration(:,i) = idx;
    [randcorr_cog(i,1) , randp_cog(i,1)] = corr(mean(DMNarray(:,idx))',allfactors(idx,16));
    [randcorr_pfactor(i,1) , randp_pfactor(i,1)] = corr(mean(DMNarray(:,idx))',allfactors(idx,30));
end

% Make figures
[a b] = sort(randcorr_cog,'descend');
[y z] = sort(randcorr_cog);
idx = randperm(10);

figure; hold on 

% Full sample
e = scatter(mean(DMNarray,1)',allfactors(:,16),0.01);
d = lsline;
xd = d.XData;
xy = d.YData; 
close all 

% subsamples + full sample 
figure; hold on
% Full sample 
d = plot([xd],[xy]);
d.Color = [0 0 0];
d.LineWidth = 4;
d.LineStyle = '--';
ylim([70 150])
xlim([.15 .50])

% Subsample
c = scatter(mean(DMNarray(:,iteration(:,b(idx(1)))))',allfactors(iteration(:,b(idx(1))),16),125,'k','filled');
c.MarkerEdgeColor = [0 0 0];
c.MarkerFaceColor = [0 109 44]./255;
c.LineWidth = 1;
a = lsline;
a.Color = [0 109 44]./255;
a.LineWidth = 4;

e = scatter(mean(DMNarray(:,iteration(:,z(idx(1)))))',allfactors(iteration(:,z(idx(1))),16),125,'k','filled');
e.Marker = 'd';
e.MarkerEdgeColor = [0 0 0];
e.MarkerFaceColor = [199 233 192]./255;
e.LineWidth = 1;
f = lsline;
f.Color = [199 233 192]./255;
f.LineWidth = 4;


% Psychopathology 
% Make figures
[a b] = sort(randcorr_pfactor,'descend');
[y z] = sort(randcorr_pfactor);
idx = randperm(10);

figure; hold on 

% Full sample
e = scatter(mean(DMNarray,1)',allfactors(:,30),0.01);
d = lsline;
xd = d.XData;
xy = d.YData; 
close all 

% subsamples + full sample 
figure; hold on
% Full sample 
d = plot([xd],[xy]);
d.Color = [0 0 0];
d.LineWidth = 4;
d.LineStyle = '--';
ylim([20 70])
xlim([.15 .50])

% Subsample
c = scatter(mean(DMNarray(:,iteration(:,b(idx(1)))))',allfactors(iteration(:,b(idx(1))),30),125,'k','filled');
c.MarkerEdgeColor = [0 0 0];
c.MarkerFaceColor = [84 39 143]./255;
c.LineWidth = 1;
a = lsline;
a.Color = [84 39 143]./255;
a.LineWidth = 4;

e = scatter(mean(DMNarray(:,iteration(:,z(idx(1)))))',allfactors(iteration(:,z(idx(1))),30),125,'k','filled');
e.Marker = 'd';
e.MarkerEdgeColor = [0 0 0];
e.MarkerFaceColor = [218 218 235]./255;
e.LineWidth = 1;
f = lsline;
f.Color = [218 218 235]./255;
f.LineWidth = 4;


%% Cortical thickness 
% Vertex level 
thickness = cell2mat(MainTable_g1.VertexThickness);
thickness_disc = reshape(thickness,59412,size(MainTable_g1,1));

thickness = cell2mat(MainTable_g2.VertexThickness);
thickness_rep = reshape(thickness,59412,size(MainTable_g2,1));
clear thickness

thickness_disc = thickness_disc(:,MainTable_g1.FrameTotal >= 600 & MainTable_g1.MissingRestSubjects == 0 & MainTable_g1.Badtmaskidx == 0);
% replication set
thickness_rep = thickness_rep(:,MainTable_g2.FrameTotal >= 600 & MainTable_g2.MissingRestSubjects == 0 & MainTable_g2.Badtmaskidx == 0);
allmats_ordered = cat(2,thickness_disc,thickness_rep);

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
allfactors(Missingthickness==1,:) = [];
allmats_ordered = allmats_ordered';

% Cognition 
% Run all brain behavior correlations
factors = allfactors(:,16); % run across these factors 
FullSampleCorrs = zeros(size(allmats_ordered,2),size(factors,2));
for x = 1:size(allmats_ordered,2)
    for f = 1:size(factors,2)
        FullSampleCorrs(x,f) = corr(squeeze(allmats_ordered(:,x)),factors(:,f));
    end
end

% Get max abs correlation
[a,idx] = max(abs(FullSampleCorrs));

% Set input params  
binsize = ([25 33 50 70 100 135 200 265 375 525 725 1000 1430 2000 2800 3604]);

iter = 1000; % how many times to draw with replacement at each interval 

% Get edge correlations resampling incrementally larger samples 
f = 16;
[Corrs,Pvals] = abcd_edgewise_correlation_iterative_reliability_single_factor(allmats_ordered(:,idx),allfactors(:,f),binsize,iter);
Corrs = squeeze(Corrs);
prctileplotting(Corrs,[116 196 118],'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/SV_ctIQ_CIs')


%% Clinical 
factors = allfactors(:,30); % run across these factors 
FullSampleCorrs = zeros(size(allmats_ordered,2),size(factors,2));
for x = 1:size(allmats_ordered,2)
    for f = 1:size(factors,2)
        FullSampleCorrs(x,f) = corr(squeeze(allmats_ordered(:,x)),factors(:,f));
    end
end

% Get max abs correlation
[a,idx] = max(abs(FullSampleCorrs));

% Set input params  
binsize = ([25 33 50 70 100 135 200 265 375 525 725 1000 1430 2000 2800 3604]);
iter = 1000; % how many times to draw with replacement at each interval 

% Get edge correlations resampling incrementally larger samples 
f = 30;
[Corrs,Pvals] = abcd_edgewise_correlation_iterative_reliability_single_factor(allmats_ordered(:,idx),allfactors(:,f),binsize,iter);
Corrs = squeeze(Corrs);
prctileplotting(Corrs,[158 154 200],'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/SV_ctPF_CIs')


%% PANEL D - DIVERGENT SCATTER PLOTS
allmats_ordered = allmats_ordered';
template = ft_read_cifti_mod('/data/nil-bluearc/GMT/Evan/MSC/Analysis_V2/Networks_template.dscalar.nii');
DMNverts = find(template.data(1:59412,1)==1);

% Grab DMN thickness for each subject 
DMNarray = allmats_ordered(DMNverts,:);
  
mDMN = mean(DMNarray,1)';
[r p] = corr(mDMN,allfactors(:,16));
disp(['Full sample DMN thickness Cog Abil correlation = ' num2str(r) ' : p value = ' num2str(p)]);
[r p] = corr(mDMN,allfactors(:,30));
disp(['Full sample DMN thickness p factor correlation = ' num2str(r) ' : p value = ' num2str(p)]);

% Loop thru 1000 iterations of this correlation, grabbing 50 subs each time
% Save out subject indices at each loop 
for i = 1:1000
    idx = datasample(1:size(DMNarray,2),25)';
    iteration(:,i) = idx;
    [randcorr_cog(i,1) randp_cog(i,1)] = corr(mean(DMNarray(:,idx))',allfactors(idx,16));
    [randcorr_pfactor(i,1) randp_pfactor(i,1)] = corr(mean(DMNarray(:,idx))',allfactors(idx,30));
end

[a b] = sort(randcorr_cog,'descend');
[y z] = sort(randcorr_cog);
idx = randperm(10);

figure; hold on 

% Full sample
e = scatter(mDMN,allfactors(:,16),0.01);
d = lsline;
xd = d.XData;
xy = d.YData; 
close all 

% subsamples + full sample 
figure; hold on
% Full sample 
d = plot([xd],[xy]);
d.Color = [0 0 0];
d.LineWidth = 4;
d.LineStyle = '--';
ylim([70 150])
xlim([3 3.5])

% Subsample
c = scatter(mean(DMNarray(:,iteration(:,b(idx(1)))))',allfactors(iteration(:,b(idx(1))),16),125,'k','filled');
c.MarkerEdgeColor = [0 0 0];
c.MarkerFaceColor = [0 109 44]./255;
c.LineWidth = 1;
g = lsline;
g.Color = [0 109 44]./255;
g.LineWidth = 4;

e = scatter(mean(DMNarray(:,iteration(:,z(idx(1)))))',allfactors(iteration(:,z(idx(1))),16),125,'k','filled');
e.Marker = 'd';
e.MarkerEdgeColor = [0 0 0];
e.MarkerFaceColor = [199 233 192]./255;
e.LineWidth = 1;
f = lsline;
f.Color = [199 233 192]./255;
f.LineWidth = 4;



% Psychopathology 

%
% Make figures
[a b] = sort(randcorr_pfactor,'descend');
[y z] = sort(randcorr_pfactor);

figure; hold on 

% Full sample
e = scatter(mean(DMNarray,1)',allfactors(:,30),0.01);
d = lsline;
xd = d.XData;
xy = d.YData; 
close all 

% subsamples + full sample 
figure; hold on
% Full sample 
d = plot([xd],[xy]);
d.Color = [0 0 0];
d.LineWidth = 4;
d.LineStyle = '--';
ylim([20 70])
xlim([3 3.5])

% Subsample
c = scatter(mean(DMNarray(:,iteration(:,b(idx(1)))))',allfactors(iteration(:,b(idx(1))),30),125,'k','filled');
c.MarkerEdgeColor = [0 0 0];
c.MarkerFaceColor = [84 39 143]./255;
c.LineWidth = 1;
g = lsline;
g.Color = [84 39 143]./255;
g.LineWidth = 4;

e = scatter(mean(DMNarray(:,iteration(:,z(idx(1)))))',allfactors(iteration(:,z(idx(1))),30),125,'k','filled');
e.Marker = 'd';
e.MarkerEdgeColor = [0 0 0];
e.MarkerFaceColor = [218 218 235]./255;
e.LineWidth = 1;
f = lsline;
f.Color = [218 218 235]./255;
f.LineWidth = 4;


%% Supplemental Figure 
% Goal is to show vibration effect for all Toolbox and CBCL; colored
% appropriately 

%% RSFC
% Lod resampling correlations 
load('/data/nil-bluearc/GMT/Scott/ABCD/Edgewise_func_corrs_43vars_4k_08FD.mat')


% Toolbox
cmap = [];
cmap(1,:) = [0 68 27]./255;
cmap(2,:) = [65 171 93]./255;
cmap(3,:) = [199 233 192]./255;
[x,y]=meshgrid([1:size(cmap,1)],[1:10]);
cmap = interp2(x([1,5,10],:),y([1,5,10],:),cmap,x,y);

figure; hold on 
for i = 2:11
    subplot(2,5,i-1)
    thesecorrs = mean(squeeze(Corrs(:,16,:,i)),2);
    [a,idx] = max(thesecorrs);
    Thesecorrs = squeeze(Corrs(idx,:,:,i));
    mthisbin = mean(Thesecorrs,2);
    lineProps.col = {[cmap(i-1,:) .5]}; % CBCL colors
    maxminthisbin = (max(Thesecorrs') - min(Thesecorrs')) ./ 2;
    mseb(1:length(mthisbin),mthisbin,maxminthisbin,lineProps,1);
%     plot(1:16,mthisbin,'Color',cmap(i-1,:),'LineWidth',3);
%     plot(1:16,max(Thesecorrs'),'Color',cmap(i-1,:),'LineWidth',3);
%     plot(1:16,min(Thesecorrs'),'Color',cmap(i-1,:),'LineWidth',3);
    xlim([1 16])
    ylim([-1 1])
end

% CBCL
cmap = [];
cmap(1,:) = [63 0 125]./255;
cmap(2,:) = [128 125 186]./255;
cmap(3,:) = [208 208 228]./255;
[x,y]=meshgrid([1:size(cmap,1)],[1:10]);
cmap = interp2(x([1,5,10],:),y([1,5,10],:),cmap,x,y);

figure; hold on 
for i = 12:21
    subplot(2,5,i-11)
    thesecorrs = mean(squeeze(Corrs(:,16,:,i)),2);
    [a,idx] = max(thesecorrs);
    Thesecorrs = squeeze(Corrs(idx,:,:,i));
    mthisbin = mean(Thesecorrs,2);
    lineProps.col = {[cmap(i-11,:) .5]}; % CBCL colors
    maxminthisbin = (max(Thesecorrs') - min(Thesecorrs')) ./ 2;
    mseb(1:length(mthisbin),mthisbin,maxminthisbin,lineProps,1);
%     plot(1:16,mthisbin,'Color',cmap(i-1,:),'LineWidth',3);
%     plot(1:16,max(Thesecorrs'),'Color',cmap(i-1,:),'LineWidth',3);
%     plot(1:16,min(Thesecorrs'),'Color',cmap(i-1,:),'LineWidth',3);
    xlim([1 16])
    ylim([-1 1])
end



%% Thickness 
% Lod resampling correlations 
load('/data/nil-bluearc/GMT/Scott/ABCD/VertThicknessCorrs_allfactors_43_4k.mat')


% Toolbox
cmap = [];
cmap(1,:) = [0 68 27]./255;
cmap(2,:) = [65 171 93]./255;
cmap(3,:) = [199 233 192]./255;
[x,y]=meshgrid([1:size(cmap,1)],[1:10]);
cmap = interp2(x([1,5,10],:),y([1,5,10],:),cmap,x,y);

figure; hold on 
for i = 2:11
    subplot(2,5,i-1)
    thesecorrs = mean(squeeze(Corrs(:,16,:,i)),2);
    [a,idx] = max(thesecorrs);
    Thesecorrs = squeeze(Corrs(idx,:,:,i));
    mthisbin = mean(Thesecorrs,2);
    lineProps.col = {[cmap(i-1,:) .5]}; % CBCL colors
    maxminthisbin = (max(Thesecorrs') - min(Thesecorrs')) ./ 2;
    mseb(1:length(mthisbin),mthisbin,maxminthisbin,lineProps,1);
%     plot(1:16,mthisbin,'Color',cmap(i-1,:),'LineWidth',3);
%     plot(1:16,max(Thesecorrs'),'Color',cmap(i-1,:),'LineWidth',3);
%     plot(1:16,min(Thesecorrs'),'Color',cmap(i-1,:),'LineWidth',3);
    xlim([1 16])
    ylim([-1 1])
end

% CBCL
cmap = [];
cmap(1,:) = [63 0 125]./255;
cmap(2,:) = [128 125 186]./255;
cmap(3,:) = [208 208 228]./255;
[x,y]=meshgrid([1:size(cmap,1)],[1:10]);
cmap = interp2(x([1,5,10],:),y([1,5,10],:),cmap,x,y);

figure; hold on 
for i = 12:21
    subplot(2,5,i-11)
    thesecorrs = mean(squeeze(Corrs(:,16,:,i)),2);
    [a,idx] = max(thesecorrs);
    Thesecorrs = squeeze(Corrs(idx,:,:,i));
    mthisbin = mean(Thesecorrs,2);
    lineProps.col = {[cmap(i-11,:) .5]}; % CBCL colors
    maxminthisbin = (max(Thesecorrs') - min(Thesecorrs')) ./ 2;
    mseb(1:length(mthisbin),mthisbin,maxminthisbin,lineProps,1);
%     plot(1:16,mthisbin,'Color',cmap(i-1,:),'LineWidth',3);
%     plot(1:16,max(Thesecorrs'),'Color',cmap(i-1,:),'LineWidth',3);
%     plot(1:16,min(Thesecorrs'),'Color',cmap(i-1,:),'LineWidth',3);
    xlim([1 16])
    ylim([-1 1])
end


