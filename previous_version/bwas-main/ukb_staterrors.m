%% Uk Biobank error curves (Fig 3 supplement)

load('/data/nil-bluearc/GMT/Scott/UKB/BWAS/ukb_rsfcmatrices_32572.mat');
load('/data/nil-bluearc/GMT/Scott/UKB/BWAS/behavioraltable_FIQ_32572.mat');

%% Load full sample correlations and pvalues for all brain-phenotype

for x = 1:size(mats,2)
    [FullSampleCorrs(x,1) , FullSamplePvals(x,1)] = corr(mats(:,x),subjecttable.FIQ);
end

%% Bootstrap subsampling procedure on all edges across all phenotypes  
iter = 100;
binsize = ([25 33 50 70 100 135 200 265 375 525 725 1000 1430 2000 2800 4000 10000 20000 30000]);
% type1 = nan(length(binsize),7);
% type2error = nan(length(binsize),7);
% typemerror = nan(length(binsize),19);
% typeserror = nan(size(FullSampleCorrs,1),19);
% typemerror_pval = cell(1,1);
% typeserror_pval = cell(1,1);

% for f = 1:size(allfactors,2)
    % Make sampling variability curves edge-wise 
    % edge 9669 is max rsfc/iq
    % edge  6123 is max rsfc/p factor 
%     disp(['On phenotype ... ' num2str(f)])
    tic;
    [Corrs,Pvals] = abcd_edgewise_correlation_iterative_reliability_single_factor(mats,subjecttable.FIQ,binsize,iter);
    toc;
    % Calculate error rates 
    [type1,type2,typem,types,types_pvals] = abcd_statisticalerrors(Corrs,Pvals,FullSampleCorrs,FullSamplePvals,iter);

% end
% Save 
save('/data/nil-bluearc/GMT/Scott/UKB/BWAS/type1error.mat','type1')
save('/data/nil-bluearc/GMT/Scott/UKB/BWAS/type2error.mat','type2')
save('/data/nil-bluearc/GMT/Scott/UKB/BWAS/typemerror.mat','typem')
save('/data/nil-bluearc/GMT/Scott/UKB/BWAS/typeserror.mat','types')
save('/data/nil-bluearc/GMT/Scott/UKB/BWAS/typeserror_pval.mat','types_pvals')

%%
mtype1 = type1;
mtype2 = type2;
typemerror_pval = typem;

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
    h = plot(1:size(mtype1,1),mtype1(:,i),'Color',cmap(i,:),'LineWidth',2);
end

xlim([1 size(mtype1,1)])
xticks([1:1:size(mtype1,1)]) 
xticklabels({'25'; '33' ; '50' ; '70' ; '100' ; '135' ; '200' ; '265' ; '375' ; '519' ; '725' ; '1000' ; '1430' ; '2000' ; '2800'; '4000' ; '10000' ; '20000' ; '30000'});
xtickangle(45)
xlabel('Number of Subjects','FontWeight','bold','FontSize',14)
ylabel('Type I Error','FontWeight','bold','FontSize',14)
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
ylim([0 100])
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/Type1error_UKB','tiffn')


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
    h = plot(1:size(mtype2,1),mtype2(:,i),'Color',cmap(i,:),'LineWidth',2);
end

xlim([1 size(mtype2,1)])
xticks([1:1:size(mtype2,1)]) 
xticklabels({'25'; '33' ; '50' ; '70' ; '100' ; '135' ; '200' ; '265' ; '375' ; '519' ; '725' ; '1000' ; '1430' ; '2000' ; '2800'; '4000' ; '10000' ; '20000' ; '30000'});
xtickangle(45)
xlabel('Number of Subjects','FontWeight','bold','FontSize',14)
ylabel('Type II Error','FontWeight','bold','FontSize',14)
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
ylim([0 100])
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/Type2error_UKB','tiffn')


% Type M
% Type M plotting
inflations = [4 9 18];
figure; hold on 
% plot the uncorrected line
%h=plot([1:16],typem_uncorrected_both','-.k*','LineWidth',2);

mtypem = [];    
for f = 1:size(typemerror_pval,1)
    AllInflatedCorrs = typemerror_pval{f};
    for x = 1:length(inflations)
       this(:,x) = nanmean(AllInflatedCorrs(:,:,inflations(x)),1)';
    end
    mtypem(:,:,f) = this;
end

mtypem(isnan(mtypem)) = 100;
    
thr = size(mtypem,1) + 5;
cmap = [];
cmap(1,:) = [237 248 251]./255;
cmap(2,:) = [140 150 198]./255;
cmap(3,:) = [110 1 107]./255;
[x,y]=meshgrid([1:3],[1:thr]);
cmap = interp2(x([1,round(thr/2),thr],:),y([1,round(thr/2),thr],:),cmap,x,y);
    
% p < 0.05
% h
% 
% % p < 0.000005
% lineProps.col = {cmap(ti,:)}; % CBCL colors
% lineProps.style = '--';
% mseb(1:size(mtypem,1),mtypem(:,1)',stypem(:,1)',lineProps,1);

figure; hold on 
for p = [1 7]
    for i = 1:length(inflations)
        if p == 1
            h = plot(1:size(mtypem,1),mtypem(:,i,p),'Color',cmap(inflations(i),:),'LineWidth',2,'LineStyle','--');
        else
            h = plot(1:size(mtypem,1),mtypem(:,i,p),'Color',cmap(inflations(i),:),'LineWidth',2);
        end
    end
end

xlim([1 size(mtypem,1)])
xticks([1:1:size(mtypem,1)]) 
xticklabels({'25'; '33' ; '50' ; '70' ; '100' ; '135' ; '200' ; '265' ; '375' ; '519' ; '725' ; '1000' ; '1430' ; '2000' ; '2800'; '4000' ; '10000' ; '20000' ; '30000'});
xtickangle(45)
xlabel('Number of Subjects','FontWeight','bold','FontSize',14)
ylabel('Type M Error','FontWeight','bold','FontSize',14)
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
ylim([0 100])
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/TypeMerror_UKB','tiffn')


% Type S
mtypes = [];
for f = 1:size(types_pvals,1)
    these = types_pvals{f};
     mtypes(:,f) = nanmean(these,1); 
end

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
% Plot uncorrected line 
plot(1:size(mtypes,1),nanmean(types,1),'-.k*','LineWidth',2);
for i = 1:7
    h = plot(1:size(mtypes,1),mtypes(:,i),'Color',cmap(i,:),'LineWidth',2);
end

xlim([1 size(mtypem,1)])
xticks([1:1:size(mtypem,1)]) 
xticklabels({'25'; '33' ; '50' ; '70' ; '100' ; '135' ; '200' ; '265' ; '375' ; '519' ; '725' ; '1000' ; '1430' ; '2000' ; '2800'; '4000' ; '10000' ; '20000' ; '30000'});
xtickangle(45)
xlabel('Number of Subjects','FontWeight','bold','FontSize',14)
ylabel('Type S Error','FontWeight','bold','FontSize',14)
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
ylim([0 100])
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/TypeSerror_UKB','tiffn')

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
    h = plot(1:size(mtype2,1),100-mtype2(:,i),'Color',cmap(i,:),'LineWidth',2);
end

xlim([1 size(mtype2,1)])
xticks([1:1:size(mtype2,1)]) 
xticklabels({'25'; '33' ; '50' ; '70' ; '100' ; '135' ; '200' ; '265' ; '375' ; '519' ; '725' ; '1000' ; '1430' ; '2000' ; '2800'; '4000' ; '10000' ; '20000' ; '30000'});
xtickangle(45)
xlabel('Number of Subjects','FontWeight','bold','FontSize',14)
ylabel('Power','FontWeight','bold','FontSize',14)
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
ylim([0 100])
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/Power_UKB','tiffn')



%% Fig Sx3E
%addpath('/data/nil-bluearc/GMT/Scott/MSC_Subcortical/Scripts/')
pthr = [.05 .01 .001 .0001 .00001 .000001 .0000001];
binsize = ([25 33 50 70 100 135 200 265 375 525 725 1000 1430 2000 2800 4000 10000 15000]);

iter = 100;

discrmats = mats(1:16286,:);
reprmats = mats(16287:end,:);
discbeh = subjecttable.FIQ(1:16286,:);
repbeh = subjecttable.FIQ(16287:end,:);

output_r_disc = nan(size(discrmats,2),iter,length(binsize));
output_p_disc = nan(size(discrmats,2),iter,length(binsize));
output_r_rep = nan(size(reprmats,2),iter,length(binsize));
output_p_rep = nan(size(reprmats,2),iter,length(binsize));

rdisc = nan(size(discrmats,2),iter,length(binsize));
pdisc = nan(size(reprmats,2),iter,length(binsize));
rrep = nan(size(discrmats,2),iter,length(binsize));
prep = nan(size(reprmats,2),iter,length(binsize));

numpools = 15;
z = parpool(numpools);
parfor b = 1:length(binsize)
    disp(['On binsize ' num2str(b)])
    bin = binsize(b);
    thisdiscr = nan(size(reprmats,2),iter); 
    thisdiscp = nan(size(reprmats,2),iter); 
    thisrepr = nan(size(reprmats,2),iter); 
    thisrepp = nan(size(reprmats,2),iter); 
    for i = 1:iter
        idxd = datasample(1:size(discrmats,1),bin)';
        idxr = datasample(1:size(reprmats,1),bin)';
        for v = 1:size(discrmats,2)
            [thisdiscr(v,i),thisdiscp(v,i)] = corr(discrmats(idxd,v),discbeh(idxd));
            [thisrepr(v,i),thisrepp(v,i)] = corr(reprmats(idxr,v),repbeh(idxr));
        end
    end
    rdisc(:,:,b) = thisdiscr;
    pdisc(:,:,b) = thisdiscp;
    rrep(:,:,b) = thisrepr;
    prep(:,:,b) = thisrepp;
end

output_r_disc= rdisc;
output_p_disc = pdisc;
output_r_rep= rrep;
output_p_rep= prep;


delete(z);

percentrep = nan(iter,length(binsize),length(pthr));
for p = 1:length(pthr)
    for i = 1:100
        for b = 1:length(binsize)
            thisrep = squeeze(output_p_rep(:,i,b));
            thisdisc = squeeze(output_p_disc(:,i,b));
            thisdiscidx = logical(thisdisc < pthr(p));
            thisrepidx = logical(thisrep < pthr(p));
            combined = thisdiscidx + thisrepidx;
            percentrep(i,b,p) = length(find(combined==2))/length(find(thisdiscidx==1));
        end
    end
end
%mean_percentrep = nanmean(percentrep,4);
%std_percentrep = nanstd(percentrep,0,4);

thr = length(pthr);
cmap = [];
cmap(3,:) = [1 100 80]./255;
cmap(2,:) = [103 169 207]./255;
cmap(1,:) = [208 209 230]./255;
[x,y]=meshgrid([1:3],[1:thr]);
cmap = interp2(x([1,round(thr/2),thr],:),y([1,round(thr/2),thr],:),cmap,x,y);

figure; hold on 
for p = 1:length(pthr)
    this = percentrep(:,:,p);
    lineProps.col = {cmap(p,:)};
    mseb(1:1:length(binsize),mean(this,1),std(this,0,1),lineProps,1);
end


xlim([1 length(binsize)+1])
xticks([1:1:length(binsize)+1])
xticklabels({'25'; '33' ; '50' ; '70' ; '100' ; '135' ; '200' ; '265' ; '375' ; '519' ; '725' ; '1000' ; '1430' ; '2000' ; '2800'; '4000' ; '10000' ; '15000' ; '30000'});
xtickangle(45)
xlabel('Sample Size','FontWeight','bold','FontSize',14)
ylabel('Percent Replication','FontWeight','bold','FontSize',14)
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
ylim([0 1])
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/UKB_Percent_rep_by_sample_size','tiffn')


