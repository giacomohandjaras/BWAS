abcd_make_big_data_table_FD_08

abcd_make_brainbehavior_variables

% generate correlations and p-values in replication sample
% ("FullSampleCorrs" from previous analysis using full sample as ground
% truth
FullSampleCorrs = nan(size(allmats_2d,2),size(allfactors,2));
FullSamplePvals = nan(size(allmats_2d,2),size(allfactors,2));

for f = 1:size(allfactors,2)
    thisbehav = allfactors(1965:3928,f);
    thismat = allmats_2d(1965:3928,:);
    thismat(isnan(thisbehav),:) = [];
    thisbehav(isnan(thisbehav)) = [];
    for r = 1:size(allmats_2d,2)   
        [FullSampleCorrs(r,f) , FullSamplePvals(r,f)] = corr(thismat(:,r),thisbehav);
    end
end

%% Bootstrap subsampling procedure on all edges across all phenotypes  
iter = 100;
binsize = [25 33 45 60 80 100 145 200 256 350 460 615 825 1100 1475 1964];
type1error = nan(length(binsize),7,size(allfactors,2));
type2error = nan(length(binsize),7,size(allfactors,2));
typeserror = nan(size(FullSampleCorrs,1),16,size(allfactors,2));
typemerror_pval = cell(size(allfactors,2),1);
typeserror_pval = cell(size(allfactors,2),1);

for f = 1:size(allfactors,2)
    % Generate corrs and pvals on discovery (in) sample
    disp(['On phenotype ... ' num2str(f)])
    [Corrs,Pvals] = abcd_edgewise_correlation_iterative_reliability_single_factor(allmats_2d(1:1964,:),allfactors(1:1964,f),binsize,iter);
    
    % Calculate error rates 
    [thist1,thist2,thistm,thists,thists_pval] = abcd_statisticalerrors(Corrs,Pvals,FullSampleCorrs(:,f),FullSamplePvals(:,f),iter);
    type1error(:,:,f) = thist1;
    type2error(:,:,f) = thist2;
    typemerror_pval{f} = thistm;
    typeserror(:,:,f) = thists;
    typeserror_pval{f} = thists_pval;
    clear Corrs Pvals 

end
% Save 
save('/data/nil-bluearc/GMT/Scott/ABCD/type1error_41vars_rsfc_pval_oos.mat','type1error')
save('/data/nil-bluearc/GMT/Scott/ABCD/type2error_41vars_rsfc_pval_oos.mat','type2error')
save('/data/nil-bluearc/GMT/Scott/ABCD/typemerror_41vars_rsfc_pvals_oos.mat','typemerror_pval')
save('/data/nil-bluearc/GMT/Scott/ABCD/typeserror_41vars_rsfc_pval_oos.mat','typeserror_pval')


%% Cortical thickness 

abcd_make_big_data_table_FD_08

abcd_make_brainbehavior_variables

thickness = cell2mat(MainTable_g1.VertexThickness);
thickness_disc = reshape(thickness,59412,size(MainTable_g1,1));

thickness = cell2mat(MainTable_g2.VertexThickness);
thickness_rep = reshape(thickness,59412,size(MainTable_g2,1));
clear thickness

thickness_disc = thickness_disc(:,MainTable_g1.FrameTotal >= 600 & MainTable_g1.MissingRestSubjects == 0 & MainTable_g1.Badtmaskidx == 0);
% replication set
thickness_rep = thickness_rep(:,MainTable_g2.FrameTotal >= 600 & MainTable_g2.MissingRestSubjects == 0 & MainTable_g2.Badtmaskidx == 0);
allmats = cat(2,thickness_disc,thickness_rep);

% Index missing subject 
Missingthickness = zeros(size(allmats,2),1);
for s = 1:length(Missingthickness)
   thissub = allmats(:,s);
   if sum(isnan(thissub)) > 0
       Missingthickness(s,1) = 1;
   else
       Missingthickness(s,1) = 0;
   end
end
% Remove any missing subjects from further analysis 
allmats(:,Missingthickness == 1) = [];
allfactors(Missingthickness == 1,:) = [];
allmats = allmats';

%% Generate full sample correlations and pvalues for all brain-phenotype
FullSampleCorrs = nan(size(allmats,2),size(allfactors,2));
FullSamplePvals = nan(size(allmats,2),size(allfactors,2));

for f = 1:size(allfactors,2)
    thisbehav = allfactors(1815:3604,f);
    thismat = allmats(1815:3604,:);
    thismat(isnan(thisbehav),:) = [];
    thisbehav(isnan(thisbehav)) = [];
    for r = 1:size(allmats,2)   
        [FullSampleCorrs(r,f) , FullSamplePvals(r,f)] = corr(thismat(:,r),thisbehav);
    end
end


%% Bootstrap subsampling procedure on all edges across all phenotypes  

iter = 100;
binsize = [25 33 45 60 80 100 145 200 256 350 460 615 825 1100 1475 1814];
type1error = nan(length(binsize),7,size(allfactors,2));
type2error = nan(length(binsize),7,size(allfactors,2));
typemerror = nan(length(binsize),19,size(allfactors,2));
typeserror = nan(size(FullSampleCorrs,1),16,size(allfactors,2));
typemerror_pval = cell(size(allfactors,2),1);
typeserror_pval = cell(size(allfactors,2),1);

for f = 1:size(allfactors,2)

    disp(['On phenotype ... ' num2str(f)])
    [Corrs,Pvals] = abcd_edgewise_correlation_iterative_reliability_single_factor(allmats,allfactors(:,f),binsize,iter);
    % Calculate error rates 
    [thist1,thist2,thistm,thists,thists_pval] = abcd_statisticalerrors(Corrs,Pvals,FullSampleCorrs(:,f),FullSamplePvals(:,f),iter);
    type1error(:,:,f) = thist1;
    type2error(:,:,f) = thist2;
    typemerror_pval{f} = thistm;
    typeserror(:,:,f) = thists;
    typeserror_pval{f} = thists_pval;
    clear Corrs Pvals 

end
% Save 
save('/data/nil-bluearc/GMT/Scott/ABCD/type1error_41vars_ct_pval_oos.mat','type1error')
save('/data/nil-bluearc/GMT/Scott/ABCD/type2error_41vars_ct_pval_oos.mat','type2error')
save('/data/nil-bluearc/GMT/Scott/ABCD/typemerror_41vars_ct_pvals_oos.mat','typemerror_pval')
save('/data/nil-bluearc/GMT/Scott/ABCD/typeserror_41vars_ct_pval_oos.mat','typeserror_pval')


%% Load rsfc and ct 

% ct 
load('/data/nil-bluearc/GMT/Scott/ABCD/type1error_41vars_ct_pval_oos.mat');
load('/data/nil-bluearc/GMT/Scott/ABCD/type2error_41vars_ct_pval_oos.mat');
load('/data/nil-bluearc/GMT/Scott/ABCD/typemerror_41vars_ct_pvals_oos.mat');
load('/data/nil-bluearc/GMT/Scott/ABCD/typeserror_41vars_ct_pval_oos.mat');

type1error_ct = type1error;
type2error_ct = type2error;
typemerror_ct = typemerror_pval;
typeserror_ct = typeserror_pval;

% rsfc
load('/data/nil-bluearc/GMT/Scott/ABCD/type1error_41vars_rsfc_pval_oos.mat');
load('/data/nil-bluearc/GMT/Scott/ABCD/type2error_41vars_rsfc_pval_oos.mat');
load('/data/nil-bluearc/GMT/Scott/ABCD/typemerror_41vars_rsfc_pvals_oos.mat');
load('/data/nil-bluearc/GMT/Scott/ABCD/typeserror_41vars_rsfc_pval_oos.mat');

type1error_rsfc = type1error;
type2error_rsfc = type2error;
typemerror_rsfc = typemerror_pval;
typeserror_rsfc = typeserror_pval;

type1error = cat(3,type1error_ct,type1error_rsfc);
type2error = cat(3,type2error_ct,type2error_rsfc);
typemerror_pval = [typemerror_ct ; typemerror_rsfc];
types_pvals = [typeserror_ct ; typeserror_rsfc];

clear *rsfc *ct

%% Get mean and std across 41 vars 
mtype1 = nanmean(type1error,3);
mtype2 = nanmean(type2error,3);


% std 
stype1 = nanstd(type1error,0,3);
stype2 = nanstd(type2error,0,3);


%% Type 1 plotting - rsfc
thr = size(mtype1,2);
cmap = [];
cmap(1,:) = [1 100 80]./255;
cmap(2,:) = [103 169 207]./255;
cmap(3,:) = [208 209 230]./255;
[x,y]=meshgrid([1:3],[1:thr]);
cmap = interp2(x([1,round(thr/2),thr],:),y([1,round(thr/2),thr],:),cmap,x,y);

figure; hold on 
for i = 1:thr
    lineProps.col = {cmap(i,:)}; % CBCL colors
    mseb(1:size(mtype1,1),mtype1(:,i)',stype1(:,i)',lineProps,1);
    %h = plot(1:size(mtype1,1),mtype1(:,i),'Color',cmap(i,:),'LineWidth',2);
end

xlim([1 size(mtype1,1)])
xticks([1:1:size(mtype1,1)]) 
xticklabels({'25'; '33' ; '50' ; '70' ; '100' ; '135' ; '200' ; '265' ; '375' ; '519' ; '725' ; '1000' ; '1430' ; '2000' ; '2800'; '3928'});
xtickangle(45)
xlabel('Number of Subjects','FontWeight','bold','FontSize',14)
ylabel('Type I Error','FontWeight','bold','FontSize',14)
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
ylim([0 22])
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/Type1error_41vars_edge_both','tiffn')


% Type 2 plotting 
thr = size(mtype2,2);
cmap = [];
cmap(1,:) = [1 100 80]./255;
cmap(2,:) = [103 169 207]./255;
cmap(3,:) = [208 209 230]./255;
[x,y]=meshgrid([1:3],[1:thr]);
cmap = interp2(x([1,round(thr/2),thr],:),y([1,round(thr/2),thr],:),cmap,x,y);

figure; hold on 
for i = 1:thr
    lineProps.col = {cmap(i,:)}; % CBCL colors
    mseb(1:size(mtype2,1),mtype2(:,i)',stype2(:,i)',lineProps,1);
    %h = plot(1:size(mtype2,1),mtype2(:,i),'Color',cmap(i,:),'LineWidth',2);
end

xlim([1 size(mtype2,1)])
xticks([1:1:size(mtype2,1)]) 
xticklabels({'25'; '33' ; '50' ; '70' ; '100' ; '135' ; '200' ; '265' ; '375' ; '519' ; '725' ; '1000' ; '1430' ; '2000' ; '2800'; '3928'});
xtickangle(45)
xlabel('Number of Subjects','FontWeight','bold','FontSize',14)
ylabel('Type II Error','FontWeight','bold','FontSize',14)
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
ylim([0 100])
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/Type2error_41vars_edge_both','tiffn')


% Type M
% Type M plotting
inflations = [4 9 16];
figure; hold on 
for i = 1:length(inflations)
    mtypem = [];
    stypem = [];
    for f = 1:size(typemerror_pval,1)
        for x = 1:size(typemerror_pval{f},1)
           AllInflatedCorrs = typemerror_pval{f};
           mtypem(:,:,x,f) = squeeze(nanmean(AllInflatedCorrs{x},1)); 
        end
    end
    stypem = nanstd(mtypem,0,4); stypem(isnan(stypem)) = 0;
    mtypem = nanmean(mtypem,4);
    
    ti = inflations(i);
    mtypem = squeeze(mtypem(:,ti,:));
    stypem = squeeze(stypem(:,ti,:));
    stypem(stypem==0) = .05;
    mtypem(isnan(mtypem)) = 100;
    
    thr = size(mtypem,1) + 5;
    cmap = [];
    cmap(1,:) = [237 248 251]./255;
    cmap(2,:) = [140 150 198]./255;
    cmap(3,:) = [110 1 107]./255;
    [x,y]=meshgrid([1:3],[1:thr]);
    cmap = interp2(x([1,round(thr/2),thr],:),y([1,round(thr/2),thr],:),cmap,x,y);
    
    % p < 0.05
    lineProps.col = {cmap(ti,:)}; % CBCL colors
    lineProps.style = '-';
    mseb(1:size(mtypem,1),mtypem(:,7)',stypem(:,7)',lineProps,1);
    
    % p < 0.000005
    lineProps.col = {cmap(ti,:)}; % CBCL colors
    lineProps.style = '--';
    mseb(1:size(mtypem,1),mtypem(:,1)',stypem(:,1)',lineProps,1);


end
xlim([1 16])
xticks([1:1:size(mtypem,1)]) 
xticklabels({'25'; '33' ; '50' ; '70' ; '100' ; '135' ; '200' ; '265' ; '375' ; '519' ; '725' ; '1000' ; '1430' ; '2000' ; '2800'; '3928'});
xtickangle(45)
xlabel('Number of Subjects','FontWeight','bold','FontSize',14)
ylabel('Inflation Rate','FontWeight','bold','FontSize',14)
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
box('off')
ylim([0 100])
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/Typemerror_41vars_edge_both_pvals','tiffn')

%%save colorbar 
figure; scatter(mtypem(1,:),100-mtypem(1,:),[],cmap);colormap(cmap);colorbar
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/Fig2ccolorbar_41vars_both','tiffn')


% Type S
mtypes = []; stypes = [];
for f = 1:size(types_pvals,1)
    these = types_pvals{f};
    for i = 1:size(these,1) 
       mtypes(:,i,f) = nanmean(these{i},1); 
       stypes(:,i,f) = nanstd(these{i},0,1);
    end
end
mtypes = nanmean(mtypes,3)';
stypes = nanmean(stypes,3)';
% Type S plotting 
thr = size(mtypes,1);
cmap = [];
cmap(1,:) = [1 100 80]./255;
cmap(2,:) = [103 169 207]./255;
cmap(3,:) = [208 209 230]./255;
[x,y]=meshgrid([1:3],[1:thr]);
cmap = interp2(x([1,round(thr/2),thr],:),y([1,round(thr/2),thr],:),cmap,x,y);
cmap(:,4) = .5;

figure; hold on 
for i = 1:thr
    lineProps.col = {cmap(i,:)}; % CBCL colors
    lineProps.style = '-';
    mseb(1:size(mtypes,2),mtypes(i,:),stypes(i,:),lineProps,1);
end

xlim([1 16])
xticks(1:1:size(mtypes,2)) 
xticklabels({'25'; '33' ; '50' ; '70' ; '100' ; '135' ; '200' ; '265' ; '375' ; '519' ; '725' ; '1000' ; '1430' ; '2000' ; '2800'; '3928'});
xtickangle(45)
xlabel('Number of Subjects','FontWeight','bold','FontSize',14)
ylabel('Sign Flip Rate','FontWeight','bold','FontSize',14)
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
box('off')
ylim([0 3])    
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/Typeserror_41vars_edge_both_pvals','tiffn')

% Power 
% Type 2 plotting 
thr = size(mtype2,2);
cmap = [];
cmap(1,:) = [1 100 80]./255;
cmap(2,:) = [103 169 207]./255;
cmap(3,:) = [208 209 230]./255;
[x,y]=meshgrid([1:3],[1:thr]);
cmap = interp2(x([1,round(thr/2),thr],:),y([1,round(thr/2),thr],:),cmap,x,y);

figure; hold on 
for i = 1:thr
    lineProps.col = {cmap(i,:)}; % CBCL colors
    mseb(1:size(mtype2,1),100-mtype2(:,i)',stype2(:,i)',lineProps,1);
    %h = plot(1:size(mtype2,1),mtype2(:,i),'Color',cmap(i,:),'LineWidth',2);
end

xlim([1 size(mtype2,1)])
xticks([1:1:size(mtype2,1)]) 
xticklabels({'25'; '33' ; '50' ; '70' ; '100' ; '135' ; '200' ; '265' ; '375' ; '519' ; '725' ; '1000' ; '1430' ; '2000' ; '2800'; '3928'});
xtickangle(45)
xlabel('Number of Subjects','FontWeight','bold','FontSize',14)
ylabel('Power','FontWeight','bold','FontSize',14)
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
ylim([0 100])
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/Power_41vars_edge_both','tiffn')

