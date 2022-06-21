%% Make figures from abcd_perc_rep.m

% Load cortical thickness 
load('Fig3Edata_replicability_corticalthickness.mat')
mean_percentrep = mean_percentrep.*100; % from decimal to percentage 
% load rsfc
load('Fig3Edata_replicability_rsfc.mat')

% combine across bootstrap samples 
all = cat(2,mean_percentrep,mean_percrep);

binsize = [25 33 45 60 80 100 145 200 256 350 460 615 825 1100 1475 1964];
pthr = [.05 .01 .001 .0001 .00001 .000001 .0000001];
thr = length(pthr);
cmap = [];
cmap(3,:) = [1 100 80]./255;
cmap(2,:) = [103 169 207]./255;
cmap(1,:) = [208 209 230]./255;
[x,y]=meshgrid([1:3],[1:thr]);
cmap = interp2(x([1,round(thr/2),thr],:),y([1,round(thr/2),thr],:),cmap,x,y);

figure; hold on 
for p = 1:length(pthr)
    this = mean_percrep(:,:,p);
    lineProps.col = {cmap(p,:)};
    mseb(1:length(binsize),nanmean(this,2)',nanstd(this,0,2)',lineProps,1);
end


xlim([1 length(binsize) + 2])
xticks([1:1:length(binsize)])
xticklabels({'25'; '33' ; '45' ; '60' ; '80' ; '100' ; '145' ; '200' ; '256' ; '350' ; '460' ; '615' ; '825' ; '1100' ; '1475'; '1964'});
xtickangle(45)
xlabel('Sample Size','FontWeight','bold','FontSize',14)
ylabel('Percent Replication','FontWeight','bold','FontSize',14)
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
ylim([0 100])
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/Percent_rep_by_sample_size_ctrsfc','tiffn');

% %% projection of replicabilty to larger samples (dotted lines in 3E)
% load('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/ABCD_pitt_washu/Scripts/Figure2/allpower.mat');
% power = allpower.ave;
% standdev = fliplr(allpower.std);
% 
% powersq = 100.*(fliplr(power.*power));
% 
% meanall = squeeze(nanmean(all,2));
% stdall = squeeze(nanstd(all,0,2));
% meanall(end+1:end+2,:) = powersq(end-1:end,:); % append 2 larger sampling bins
% stdall(end+1:end+2,:) = standdev(end-1:end,:);
% 
% figure; hold on 
% for p = 1:length(pthr)
%     lineProps.col = {cmap(p,:)};
%     lineProps.style = '-';
%     mseb(1:length(binsize),meanall(1:16,p)',stdall(1:16,p)',lineProps,1);
%     lineProps.style = '--';
%     mseb(16:18,meanall(16:18,p)',stdall(16:18,p)',lineProps,1);
% end
% 
% xlim([1 length(binsize) + 2])
% xticks([1:1:length(binsize) + 2])
% xticklabels({'25'; '33' ; '45' ; '60' ; '80' ; '100' ; '145' ; '200' ; '256' ; '350' ; '460' ; '615' ; '825' ; '1100' ; '1475'; '1964'; '2800' ; '3928'});
% xtickangle(45)
% xlabel('Sample Size','FontWeight','bold','FontSize',14)
% ylabel('Percent Replication','FontWeight','bold','FontSize',14)
% set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
% ylim([0 60])
% saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/Percent_rep_by_sample_size_ctrsfc_withinterp','tiffn');


%% Effect size similarity


