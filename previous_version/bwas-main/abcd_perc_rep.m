%% ABCD % replication varying p-value in disc and rep set 2 x 2 table 

abcd_make_big_data_table_FD_08

abcd_make_brainbehavior_variables

%% run correlations on disc and rep set 
rr = []; pd = []; rr = []; pr = [];
for f = 1:size(allfactors,2)
    for i = 1:size(allmats_2d,2)
        [rd(i,f),pd(i,f)] = corr(allmats_2d(1:1964,i),allfactors(1:1964,f));
        [rr(i,f),pr(i,f)] = corr(allmats_2d(1965:size(allmats_2d,1),i),allfactors(1965:size(allmats_2d,1),f));
    end
end
% Vectorize
pd(:,[6 18 19 37:41]) = [];
pr(:,[6 18 19 37:41]) = [];
pd = pd(:);
pr = pr(:);

%% threshold disc and rep by pavlue 

thr = [.05 .01 .001 .0001 .00001 .000001 .0000001];

for d = 1:length(thr)
    for r = 1:length(thr)
        dthr = thr(d);
        rthr = thr(r);
        thisdisc = logical(pd < dthr);
        thisrep = logical(pr < rthr);
        replications = thisdisc + thisrep; % 2 = replication 
        % find total number of significant r's in discovery set (ie. where
        % thisdisc = 1 
        totalsigdisc = length(find(thisdisc==1));
        % Determine % replication (2 in replications / totalsigndisc
        alloutput(r,d) = 100*(length(find(replications==2)) / totalsigdisc);
    end
end


%% At each sample size across 100 iterations
% RSFC
addpath('/data/nil-bluearc/GMT/Scott/MSC_Subcortical/Scripts/')
pthr = [.05 .01 .001 .0001 .00001 .000001 .0000001];
binsize = [25 33 45 60 80 100 145 200 256 350 460 615 825 1100 1475 1964];
iter = 100;

numpools = 14;
z = parpool(numpools);

discrmats = allmats_2d(1:1964,:);
reprmats = allmats_2d(1965:3928,:);
discbeh = allfactors(1:1964,:);
repbeh = allfactors(1965:3928,:);

output_r_disc = nan(size(discrmats,2),length(binsize),iter,41);
output_p_disc = nan(size(discrmats,2),length(binsize),iter,41);
output_r_rep = nan(size(reprmats,2),length(binsize),iter,41);
output_p_rep = nan(size(reprmats,2),length(binsize),iter,41);

parfor f = 1:size(allfactors,2)
    rdisc = nan(size(discrmats,2),length(binsize),iter);
    pdisc = nan(size(reprmats,2),length(binsize),iter);
    rrep = nan(size(discrmats,2),length(binsize),iter);
    prep = nan(size(reprmats,2),length(binsize),iter);
    for b = 1:length(binsize)
        disp(['On factor' num2str(f) ' binsize ' num2str(b)])
        bin = binsize(b);
        for i = 1:iter
            idxd = datasample(1:size(discrmats,1),bin)';
            idxr = datasample(1:size(reprmats,1),bin)';
            for v = 1:size(discrmats,2)
                [rdisc(v,b,i),pdisc(v,b,i)] = corr(discrmats(idxd,v),discbeh(idxd,f));
                [rrep(v,b,i),prep(v,b,i)] = corr(reprmats(idxr,v),repbeh(idxr,f));
            end
        end

    end

    output_r_disc(:,:,:,f) = rdisc;
    output_p_disc(:,:,:,f) = pdisc;
    output_r_rep(:,:,:,f) = rrep;
    output_p_rep(:,:,:,f) = prep;
end

delete(z);

percentrep = nan(length(binsize),iter,length(pthr),size(allfactors,2));
for p = 1:length(pthr)
    for f = 1:41    
        for i = 1:100
            for b = 1:16
                thisrep = squeeze(output_p_rep(:,b,i,f));
                thisdisc = squeeze(output_p_disc(:,b,i,f));
                thisdiscidx = logical(thisdisc < pthr(p));
                thisrepidx = logical(thisrep < pthr(p));
                combined = thisdiscidx + thisrepidx;
                percentrep(b,i,p,f) = length(find(combined==2))/length(find(thisdiscidx==1));
            end
        end
    end
end
mean_percentrep = nanmean(percentrep,4);
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
    this = mean_percentrep(:,:,p);
    lineProps.col = {cmap(p,:)};
    mseb(1:length(binsize),mean(this,2)',std(this,0,2)',lineProps,1);
end


xlim([1 length(binsize)])
xticks([1:1:length(binsize)])
xticklabels({'25'; '33' ; '45' ; '60' ; '80' ; '100' ; '145' ; '200' ; '256' ; '350' ; '460' ; '615' ; '825' ; '1100' ; '1475'; '1964'});
xtickangle(45)
xlabel('Sample Size','FontWeight','bold','FontSize',14)
ylabel('Percent Replication','FontWeight','bold','FontSize',14)
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
ylim([0 1])
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/Percent_rep_by_sample_size','tiffn')

% compare top 1% 
percentreduction = nan(length(binsize),iter,size(allfactors,2));
magnitudereduction = nan(length(binsize),iter,size(allfactors,2));
discmean = nan(length(binsize),iter,size(allfactors,2));
repmean = nan(length(binsize),iter,size(allfactors,2));
for f = 1:41    
    for i = 1:100
        for b = 1:16
            thisdisc = squeeze(output_r_disc(:,b,i,f));
            [a j] = sort(thisdisc,'descend');
            thisdisc = a(1:774);
            thisrep = squeeze(output_r_rep(j(1:774),b,i,f));
            discmean(b,i,f) = nanmean(thisdisc);
            repmean(b,i,f) = nanmean(thisrep);
            percentreduction(b,i,f) = nanmean(100*((thisdisc - thisrep)./thisdisc));
            magnitudereduction(b,i,f) = nanmean(thisrep-thisdisc);
        end
    end
end
mean_percentreduction = mean(nanmean(percentreduction,2),2);
mean_magnitudereduction = mean(nanmean(magnitudereduction,2),2);
mean_discmean = mean(nanmean(discmean,3),2);
mean_repmean = mean(nanmean(repmean,3),2);



%% Cortical thickness

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


pthr = [.05 .01 .001 .0001 .00001 .000001 .0000001];
binsize = [25 33 45 60 80 100 145 200 256 350 460 615 825 1100 1475 1964];
iter = 100;

numpools = 14;
z = parpool(numpools);

discrmats = allmats(1:1814,:);
reprmats = allmats(1815:3604,:);
discbeh = allfactors(1:1814,:);
repbeh = allfactors(1815:3604,:);

output_r_disc = nan(size(discrmats,2),length(binsize),iter,41);
output_p_disc = nan(size(discrmats,2),length(binsize),iter,41);
output_r_rep = nan(size(reprmats,2),length(binsize),iter,41);
output_p_rep = nan(size(reprmats,2),length(binsize),iter,41);

parfor f = 1:size(allfactors,2)
    rdisc = nan(size(discrmats,2),length(binsize),iter);
    pdisc = nan(size(reprmats,2),length(binsize),iter);
    rrep = nan(size(discrmats,2),length(binsize),iter);
    prep = nan(size(reprmats,2),length(binsize),iter);
    for b = 1:length(binsize)
        disp(['On factor' num2str(f) ' binsize ' num2str(b)])
        bin = binsize(b);
        for i = 1:iter
            idxd = datasample(1:size(discrmats,1),bin)';
            idxr = datasample(1:size(reprmats,1),bin)';
            for v = 1:size(discrmats,2)
                [rdisc(v,b,i),pdisc(v,b,i)] = corr(discrmats(idxd,v),discbeh(idxd,f));
                [rrep(v,b,i),prep(v,b,i)] = corr(reprmats(idxr,v),repbeh(idxr,f));
            end
        end

    end

    output_r_disc(:,:,:,f) = rdisc;
    output_p_disc(:,:,:,f) = pdisc;
    output_r_rep(:,:,:,f) = rrep;
    output_p_rep(:,:,:,f) = prep;
end

delete(z);

percentrep = nan(length(binsize),iter,length(pthr),size(allfactors,2));
for p = 1:length(pthr)
    for f = 1:41    
        for i = 1:100
            for b = 1:16
                thisrep = squeeze(output_p_rep(:,b,i,f));
                thisdisc = squeeze(output_p_disc(:,b,i,f));
                thisdiscidx = logical(thisdisc < pthr(p));
                thisrepidx = logical(thisrep < pthr(p));
                combined = thisdiscidx + thisrepidx;
                percentrep(b,i,p,f) = length(find(combined==2))/length(find(thisdiscidx==1));
            end
        end
    end
end
mean_percentrep = nanmean(percentrep,4);
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
    this = mean_percentrep(:,:,p);
    lineProps.col = {cmap(p,:)};
    mseb(1:length(binsize),mean(this,2)',std(this,0,2)',lineProps,1);
end


xlim([1 length(binsize)])
xticks([1:1:length(binsize)])
xticklabels({'25'; '33' ; '45' ; '60' ; '80' ; '100' ; '145' ; '200' ; '256' ; '350' ; '460' ; '615' ; '825' ; '1100' ; '1475'; '1814'});
xtickangle(45)
xlabel('Sample Size','FontWeight','bold','FontSize',14)
ylabel('Percent Replication','FontWeight','bold','FontSize',14)
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
ylim([0 .3])
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/Percent_rep_by_sample_size_ct','tiffn')

% compare top 1% 
percentreduction = nan(length(binsize),iter,size(allfactors,2));
magnitudereduction = nan(length(binsize),iter,size(allfactors,2));
discmean = nan(length(binsize),iter,size(allfactors,2));
repmean = nan(length(binsize),iter,size(allfactors,2));
for f = 1:41    
    for i = 1:100
        for b = 1:16
            thisdisc = squeeze(output_r_disc(:,b,i,f));
            [a j] = sort(thisdisc,'descend');
            thisdisc = a(1:774);
            thisrep = squeeze(output_r_rep(j(1:553),b,i,f));
            discmean(b,i,f) = nanmean(thisdisc);
            repmean(b,i,f) = nanmean(thisrep);
            percentreduction(b,i,f) = nanmean(100*((thisdisc - thisrep)./thisdisc));
            magnitudereduction(b,i,f) = nanmean(thisrep-thisdisc);
        end
    end
end
mean_percentreduction = mean(nanmean(percentreduction,2),2);
mean_magnitudereduction = mean(nanmean(magnitudereduction,2),2);
mean_discmean = mean(nanmean(discmean,3),2);
mean_repmean = mean(nanmean(repmean,3),2);

% Save to structure
effectsizesimilarity_ct.mean_percentreduction = mean_percentreduction; 
effectsizesimilarity_ct.mean_magnitudereduction = mean_magnitudereduction; 
effectsizesimilarity_ct.mean_discmean = mean_discmean; 
effectsizesimilarity_ct.mean_repmean = mean_repmean; 

save('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/effectsizesimilarity_ct.mat','effectsizesimilarity_ct');