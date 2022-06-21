% Grab top 5% of discovery correlations ... determine % replication 


%% CCA rsfc 
load('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/CCA/insampleweightsall_20percentcomps.mat');
load('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/CCA/outofsampleweightsall_20percentcomps.mat');

discovery = insampleweightsall.tlbx;
replication = outofsampleweightsall.tlbx;

discovery = cat(3,discovery,insampleweightsall.cbcl);
replication = cat(3,replication,outofsampleweightsall.cbcl);

% CCA ct 
load('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/CCA/insampleweightsall_20percentcomp_thickness.mat');
load('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/CCA/outofsampleweightsall_20percentcomp_thickness.mat');

discovery = cat(3,discovery,insampleweightsall.tlbx);
replication = cat(3,replication,outofsampleweightsall.tlbx);
discovery = cat(3,discovery,insampleweightsall.cbcl);
replication = cat(3,replication,outofsampleweightsall.cbcl);



%% svr 
rsfc = readtable('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/PCA5update20200819.rsfc.csv' , 'Delimiter',',');
ct = readtable('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/PCA5update20200819.vertex.csv','Delimiter',',');
rsfc.traincor(3200) = .5; rsfc.outofsamplecor(3200) = .35;

%ct tlbx
discovery = cat(3,discovery,reshape(ct.traincor(1:1600),100,16)');
replication = cat(3,replication,reshape(ct.outofsamplecor(1:1600),100,16)');

ct = readtable('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/PCA5_20200715vertex.csv','Delimiter',',');
% ct cbcl 
discovery = cat(3,discovery,reshape(ct.traincor(1:1600),100,16)');
replication = cat(3,replication,reshape(ct.outofsamplecor(1:1600),100,16)');

%rsfc tlbx 
discovery = cat(3,discovery,reshape(rsfc.traincor(1601:3200),100,16)');
replication = cat(3,replication,reshape(rsfc.outofsamplecor(1601:3200),100,16)');

%rsfc cbcl 
discovery = cat(3,discovery,reshape(rsfc.traincor(1:1600),100,16)');
replication = cat(3,replication,reshape(rsfc.outofsamplecor(1:1600),100,16)');

percrep = [];
for i = 1:8
   for s = 1:size(discovery,1)
       this = squeeze(discovery(s,:,i));
       [a,b] = sort(this,'descend');
       rep = replication(s,b(1:5),i);
       percrep(s,i) = 20*(length(rep(rep>0.1)));
   end
   percrep(:,i) = sort(percrep(:,i),'ascend')
end
percrepmean = mean(percrep,2)

%% Continuous rep - disc (largest sample) cca vs svr
figure; 
hold on 
subplot(1,3,1)

fs_disc_cca = squeeze(discovery(16,:,1:4));
fs_rep_cca = squeeze(replication(16,:,1:4));
h=histogram(fs_rep_cca - fs_disc_cca,30); hold on 
h.FaceColor = [.8 0 0];

fs_disc_svr = squeeze(discovery(16,:,5:8));
fs_rep_svr = squeeze(replication(16,:,5:8));
j=histogram(fs_rep_svr - fs_disc_svr,30);
j.FaceColor = [.8 .4 0];

% imaging modality
subplot(1,3,2)
hold on
fs_disc_ct = squeeze(discovery(16,:,3:6));
fs_rep_ct = squeeze(replication(16,:,3:6));
h=histogram(fs_rep_ct - fs_disc_ct,30); hold on 
h.FaceColor = [235 218 32]./255;

fs_disc_rsfc = squeeze(discovery(16,:,[1 2 7 8]));
fs_rep_rsfc = squeeze(replication(16,:,[1 2 7 8]));
j=histogram(fs_rep_rsfc - fs_disc_rsfc,30);
j.FaceColor = [107 174 214]./255;

% phenotype
subplot(1,3,3)
hold on
fs_disc_cbcl = squeeze(discovery(16,:,[2 4 6 8]));
fs_rep_cbcl = squeeze(replication(16,:,[2 4 6 8]));
h=histogram(fs_rep_cbcl - fs_disc_cbcl,30); hold on 
h.FaceColor = [158 154 200]./255;

fs_disc_tlbx = squeeze(discovery(16,:,[1 3 5 7]));
fs_rep_tlbx = squeeze(replication(16,:,[1 3 5 7]));
j=histogram(fs_rep_tlbx - fs_disc_tlbx,30);
j.FaceColor = [116 196 118]./255;

%% plot percrep 
cmap(1,:) = [0 109 44]; % cca rsfc tlbx
cmap(2,:) = [84 39 143]; % cca rsfc cbcl
cmap(3,:) = [161 217 155]; % cca ct tlbx
cmap(4,:) = [188 189 220]; % cca ct cbcl 
cmap(5,:) = [199 233 192]; % svr ct tlbx
cmap(6,:) = [218 218 235]; % svr ct cbcl 
cmap(7,:) = [65 171 93]; % svr rsfc tlbx
cmap(8,:) = [128 125 186]; % svr rsfc cbcl

figure; hold on 
for i = 1:8
    plot(1:16,percrep(:,i),'Color',cmap(i,:)./255,'LineWidth',2.5);
end

%% continuous difference - cca vs svr 

rd_diff = replication - discovery;

% CCA 
rd_diff_cca = mean(rd_diff(:,:,1:4),3);
rd_diff_svr = mean(rd_diff(:,:,5:8),3);

% plot 
lineProps.col = {[.8 0 0]};
figure;
subplot(1,3,1)
hold on 
mseb([1:16],mean(rd_diff_cca,2)',std(rd_diff_cca,0,2)',lineProps,1);
lineProps.col = {[.8 .4 0]};
mseb([1:16],mean(rd_diff_svr,2)',std(rd_diff_svr,0,2)',lineProps,1);

xlim([1 16])
xticks([1:16])
xticklabels({'25'; '35' ; '45' ; '60' ; '80' ; '100' ; '145' ; '200' ; '256' ; '350' ; '460' ; '615' ; '825' ; '1100' ; '1475'; '1964'});
xtickangle(45)
xlabel('Number of Subjects','FontWeight','bold','FontSize',14)
ylabel('Replication r - Discovery r','FontWeight','bold','FontSize',14)
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
box('off')
ylim([-1 0])


%% continuous difference - ct vs rsfc 

rd_diff = replication - discovery;

% ct & rsfc
rd_diff_ct = mean(rd_diff(:,:,[3:6]),3);
rd_diff_rsfc = mean(rd_diff(:,:,[1 2 7 8]),3);

% plot 
lineProps.col = {[.8 0 0]};
figure;
subplot(1,3,2)
hold on 
mseb([1:16],mean(rd_diff_ct,2)',std(rd_diff_ct,0,2)',lineProps,1);
lineProps.col = {[.8 .4 0]};
mseb([1:16],mean(rd_diff_rsfc,2)',std(rd_diff_rsfc,0,2)',lineProps,1);

xlim([1 16])
xticks([1:16])
xticklabels({'25'; '35' ; '45' ; '60' ; '80' ; '100' ; '145' ; '200' ; '256' ; '350' ; '460' ; '615' ; '825' ; '1100' ; '1475'; '1964'});
xtickangle(45)
xlabel('Number of Subjects','FontWeight','bold','FontSize',14)
ylabel('Replication r - Discovery r','FontWeight','bold','FontSize',14)
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
box('off')
ylim([-1 0])

%% continuous difference - cbcl vs tlbx 

rd_diff = replication - discovery;

% CCA 
rd_diff_cca = mean(rd_diff(:,:,1:4),3);
rd_diff_svr = mean(rd_diff(:,:,5:8),3);

% plot 
lineProps.col = {[.8 0 0]};
figure;
subplot(1,3,3)
hold on 
mseb([1:16],mean(rd_diff_cca,2)',std(rd_diff_cca,0,2)',lineProps,1);
lineProps.col = {[.8 .4 0]};
mseb([1:16],mean(rd_diff_svr,2)',std(rd_diff_svr,0,2)',lineProps,1);

xlim([1 16])
xticks([1:16])
xticklabels({'25'; '35' ; '45' ; '60' ; '80' ; '100' ; '145' ; '200' ; '256' ; '350' ; '460' ; '615' ; '825' ; '1100' ; '1475'; '1964'});
xtickangle(45)
xlabel('Number of Subjects','FontWeight','bold','FontSize',14)
ylabel('Replication r - Discovery r','FontWeight','bold','FontSize',14)
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
box('off')
ylim([-1 0])
