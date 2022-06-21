abcd_make_big_data_table_FD_08

abcd_make_brainbehavior_variables


%% CCA computing svd at each sample size iteration 
%allmats_2d = allmats_2d';

% Set input params 
%numcomp = 1000; % number of components 
insamp_brain = allmats_2d(1:1964,:); 
outofsamp_brain = allmats_2d(1965:end,:); 
%sd = 16; % sampling density (ie how many intervals on log scale do you want
%binsize = round(logspace(log10(25),log10(size(insamp_brain,2)),sd));
binsize = [25 33 45 60 80 100 145 200 256 350 460 615 825 1100 1475 1964];
iter = 100;

% Toolbox
insamp_behavior = allfactors(1:size(rmats_disc,3),[7:13]);
outofsamp_behavior = allfactors(1+size(rmats_disc,3):size(rmats_disc,3)+size(rmats_rep,3),[7:13]);
% run loop 
tic;
[outofsampleweights_tlbx,insampleweights_tlbx] = run_svd_cca(insamp_brain,insamp_behavior,outofsamp_brain,outofsamp_behavior,binsize,iter);
toc;

% CBCL
insamp_behavior = allfactors(1:size(rmats_disc,3),[20:27]);
outofsamp_behavior = allfactors(1+size(rmats_disc,3):size(rmats_disc,3)+size(rmats_rep,3),[20:27]);
% run loop 
tic;
[outofsampleweights_cbcl,insampleweights_cbcl] = run_svd_cca(insamp_brain,insamp_behavior,outofsamp_brain,outofsamp_behavior,binsize,iter);
toc;

% Both toolbox and cbcl 
insamp_behavior = allfactors(1:size(rmats_disc,3),[2:8 12:18]);
outofsamp_behavior = allfactors(1+size(rmats_disc,3):size(rmats_disc,3)+size(rmats_rep,3),[2:8 12:18]);
% run loop 
tic;
[outofsampleweights_both,insampleweights_both] = run_svd_cca(insamp_brain,insamp_behavior,outofsamp_brain,outofsamp_behavior,binsize,iter);
toc;

% Save output 
outofsampleweightsall.tlbx = outofsampleweights_tlbx;
outofsampleweightsall.cbcl = outofsampleweights_cbcl;
%outofsampleweightsall.both = outofsampleweights_both;
insampleweightsall.tlbx = insampleweights_tlbx;
insampleweightsall.cbcl = insampleweights_cbcl;
%insampleweightsall.both = insampleweights_both;
save('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/CCA/outofsampleweightsall_20percentcomps.mat','outofsampleweightsall');
save('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/CCA/insampleweightsall_20percentcomps.mat','insampleweightsall');

% plot
figure; hold on 
% both
Thesecorrs = outofsampleweights_both;
mthisbin = mean(Thesecorrs,2);
lineProps.col = {[125/255 125/255 125/255 .2]}; % both colors
maxminthisbin = (max(Thesecorrs') - min(Thesecorrs')) ./ 2;
mseb(1:length(mthisbin),mthisbin,maxminthisbin,lineProps,1);

% toolbox 
Thesecorrs = outofsampleweights_tlbx;
mthisbin = mean(Thesecorrs,2);
lineProps.col = {[116 196 118 255]./255}; % toolbox colors
maxminthisbin = (max(Thesecorrs') - min(Thesecorrs')) ./ 2;
mseb(1:length(mthisbin),mthisbin,maxminthisbin,lineProps,1);

% cbcl
Thesecorrs = outofsampleweights_cbcl;
mthisbin = mean(Thesecorrs,2);
lineProps.col = {[158 154 200 255]./255}; % CBCL colors
maxminthisbin = (max(Thesecorrs') - min(Thesecorrs')) ./ 2;
mseb(1:length(mthisbin),mthisbin,maxminthisbin,lineProps,1);

binsize = [25 33 45 60 80 100 145 200 256 350 460 615 825 1100 1475 1964];

xlim([1 length(binsize)])
ylim([-.2 .45])
xticklabels({num2str(binsize(2)); num2str(binsize(4)); num2str(binsize(6)); num2str(binsize(8)); num2str(binsize(10)); num2str(binsize(12)); num2str(binsize(14)); num2str(binsize(16))}); 
xlabel('Number of Subjects','FontWeight','bold','FontSize',14)
ylabel('Out of Sample Correlation (r)','FontWeight','bold','FontSize',14)
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/CCA_outsampcorr_20percentcomponents','tiffn')





%% cortical thickness vertex wise 

abcd_make_big_data_table_FD_08

abcd_make_brainbehavior_variables

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

% Run first on vertex wise matrices 
allmats_2d = allmats_ordered_v';

% Set input params 
%numcomp = 20; % number of components 
insamp_brain = allmats_2d(1:1814,:); 
outofsamp_brain = allmats_2d(1815:end,:); 
%sd = 16; % sampling density (ie how many intervals on log scale do you want
%binsize = round(logspace(log10(25),log10(size(insamp_brain,2)),sd));
binsize = [25 33 45 60 80 100 145 200 256 350 460 615 825 1100 1475 1814];
iter = 100;

% covariate - FD 
%cov_insamp = MainTable_g1.FD(MainTable_g1.FrameTotal >= 600 & MainTable_g1.MissingRestSubjects == 0 & MainTable_g1.Badtmaskidx == 0);
%cov_oos = MainTable_g2.FD(MainTable_g2.FrameTotal >= 600 & MainTable_g2.MissingRestSubjects == 0 & MainTable_g2.Badtmaskidx == 0);

% Toolbox
insamp_behavior = allfactors(1:1814,[7:13]);
outofsamp_behavior = allfactors(1815:end,[7:13]);
% run loop 
tic;
[outofsampleweights_tlbx,insampleweights_tlbx] = run_svd_cca(insamp_brain,insamp_behavior,outofsamp_brain,outofsamp_behavior,binsize,iter);
toc;

% CBCL
insamp_behavior = allfactors(1:1814,[20:27]);
outofsamp_behavior = allfactors(1815:end,[20:27]);
% run loop 
tic;
[outofsampleweights_cbcl,insampleweights_cbcl] = run_svd_cca(insamp_brain,insamp_behavior,outofsamp_brain,outofsamp_behavior,binsize,iter);
toc;

% Both toolbox and cbcl 
% insamp_behavior = allfactors(1:1814,[7:13 20:27]);
% outofsamp_behavior = allfactors(1815:end,[7:13 20:27]);
% % run loop 
% tic;
% [outofsampleweights_both,insampleweights_both] = run_svd_cca(insamp_brain,insamp_behavior,outofsamp_brain,outofsamp_behavior,binsize,iter);
% toc;

% Save output 
outofsampleweightsall.tlbx = outofsampleweights_tlbx;
outofsampleweightsall.cbcl = outofsampleweights_cbcl;
%outofsampleweightsall.both = outofsampleweights_both;
insampleweightsall.tlbx = insampleweights_tlbx;
insampleweightsall.cbcl = insampleweights_cbcl;
%insampleweightsall.both = insampleweights_both;
save('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/CCA/outofsampleweightsall_20percentcomp_thickness.mat','outofsampleweightsall');
save('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/CCA/insampleweightsall_20percentcomp_thickness.mat','insampleweightsall');

% plot
figure; hold on 
% both
Thesecorrs = outofsampleweights_both;
mthisbin = mean(Thesecorrs,2);
lineProps.col = {[125/255 125/255 125/255 .2]}; % both colors
maxminthisbin = (max(Thesecorrs') - min(Thesecorrs')) ./ 2;
mseb(1:length(mthisbin),mthisbin,maxminthisbin,lineProps,1);

% toolbox 
Thesecorrs = outofsampleweights_tlbx;
mthisbin = mean(Thesecorrs,2);
lineProps.col = {[116 196 118 255]./255}; % toolbox colors
maxminthisbin = (max(Thesecorrs') - min(Thesecorrs')) ./ 2;
mseb(1:length(mthisbin),mthisbin,maxminthisbin,lineProps,1);

% cbcl
Thesecorrs = outofsampleweights_cbcl;
mthisbin = mean(Thesecorrs,2);
lineProps.col = {[158 154 200 255]./255}; % CBCL colors
maxminthisbin = (max(Thesecorrs') - min(Thesecorrs')) ./ 2;
mseb(1:length(mthisbin),mthisbin,maxminthisbin,lineProps,1);

binsize = [25 35 45 60 80 100 145 200 256 350 460 615 825 1100 1475 1964];

xlim([1 length(binsize)])
ylim([-.2 .45])
xticklabels({num2str(binsize(2)); num2str(binsize(4)); num2str(binsize(6)); num2str(binsize(8)); num2str(binsize(10)); num2str(binsize(12)); num2str(binsize(14)); num2str(binsize(16))}); 
xlabel('Number of Subjects','FontWeight','bold','FontSize',14)
ylabel('Out of Sample Correlation (r)','FontWeight','bold','FontSize',14)
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/CCA_outsampcorr_20percentcomp_thickness','tiffn')


%% Supplemental figure (6 combos: 3(20 comps 50%, all) x 2(rsfc, ct))

abcd_make_big_data_table_FD_08

abcd_make_brainbehavior_variables

allmats_2d = allmats_2d';



%% Full sample 
% insample behavior 
tlbx_insamp_behavior = allfactors(1:1964,[7:13]);
cbcl_insamp_behavior = allfactors(1:1964,[20:27]);
both_insamp_behavior = allfactors(1:1964,[7:13 20:27]);

% out of sample behavior 
tlbx_outofsamp_behavior = allfactors(1965:end,[7:13]);
cbcl_outofsamp_behavior = allfactors(1965:end,[20:27]);
both_outofsamp_behavior = allfactors(1965:end,[7:13 20:27]);

tlbx_compnum = zeros(1963,100);
cbcl_compnum = zeros(1963,100);
both_compnum = zeros(1963,100);
varexplained = zeros(1963,100);
z = parpool(18);
parfor i = 1:100
    disp(['** On iteration : ' num2str(i) ' ** '])
    
    % bootstrapped index
    idx = randperm(3928);
    insamp_brain = allmats_2d(idx(1:1964),:);

    % out of sample components 
    outofsamp_brain = allmats_2d(idx(1965:end),:); 

    % svd in-sample 
    [coeff,scores,~,~,explained,mu] = pca(insamp_brain);
    varexplained(:,i) = cumsum(explained);

    % multiple weights to oos brain data -> same coordinate space
    pc_brain_rep = (outofsamp_brain - mu)*coeff; 


    % insample behavior 
    tlbx_insamp_behavior = allfactors(idx(1:1964),[7:13]);
    cbcl_insamp_behavior = allfactors(idx(1:1964),[20:27]);
    both_insamp_behavior = allfactors(idx(1:1964),[7:13 20:27]);

    % out of sample behavior 
    tlbx_outofsamp_behavior = allfactors(idx(1965:end),[7:13]);
    cbcl_outofsamp_behavior = allfactors(idx(1965:end),[20:27]);
    both_outofsamp_behavior = allfactors(idx(1965:end),[7:13 20:27]);
    thistlbx = nan(size(scores,2),1);
    thiscbcl = nan(size(scores,2),1);
    thisboth = nan(size(scores,2),1);
    for n = 1:size(scores,2)
        % toolbox
        [A,B,~] = canoncorr(scores(:,1:n),tlbx_insamp_behavior);
        thistlbx(n,1) = corr((pc_brain_rep(:,1:n)*A(:,1)),(tlbx_outofsamp_behavior*B(:,1)));
        % toolbox
        [A,B,~] = canoncorr(scores(:,1:n),cbcl_insamp_behavior);
        thiscbcl(n,1) = corr((pc_brain_rep(:,1:n)*A(:,1)),(cbcl_outofsamp_behavior*B(:,1)));
        % toolbox
        [A,B,~] = canoncorr(scores(:,1:n),both_insamp_behavior);
        thisboth(n,1) = corr((pc_brain_rep(:,1:n)*A(:,1)),(both_outofsamp_behavior*B(:,1)));
    end
tlbx_compnum(:,i) = thistlbx;
cbcl_compnum(:,i) = thiscbcl;
both_compnum(:,i) = thisboth;
end
delete(z)

mtlbx = mean(tlbx_compnum,2);
mcbcl = mean(cbcl_compnum,2);
mboth = mean(both_compnum,2);
mve = mean(varexplained,2);
figure; hold on
h = plot(mtlbx);
h.Color = [116 198 118]./255;
h.LineWidth  = 2.5;
h = plot(mcbcl);
h.Color = [158 154 200]./255;
h.LineWidth  = 2.5;
% both
h = plot(mboth);
h.Color = [125 125 125]./255;
h.LineWidth  = 2.5;

xlim([0 1965])
ylabel('Out of sample correlation')
xlabel('Number of Components')
box('off')
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/FigSx_cca_ooscorr_numcomps_nfs_rsfc','tiffn')

% plot variance explained 
figure; 
plot(varexplained,'Color',[.2 .2 .2],'LineWidth',2);
xlabel('Number of Components')
ylabel('Cumulative Variance Explained (%)')
ylim([0 101])
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
box('off')
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/FigSx_cca_cumvar_numcomps_nfs_rsfc','tiffn')

% save output 
cca_fs_replication_numcomps.tlbxcompnum = tlbx_compnum;
cca_fs_replication_numcomps.cbclcompnum = cbcl_compnum;
cca_fs_replication_numcomps.bothcompnum = both_compnum;
cca_fs_replication_numcomps.varexplained = varexplained;
save('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/cca_fs_replication_numcomps.mat','cca_fs_replication_numcomps')



%% 300 subjects 
% initialize 300
tlbx_compnum = zeros(299,100);
cbcl_compnum = zeros(299,100);
both_compnum = zeros(299,100);
varexplained = zeros(299,100);
z = parpool(18);
parfor i = 1:100
    disp(['** On iteration : ' num2str(i) ' ** '])
    
    % bootstrapped index
    idx = randperm(3928);
    insamp_brain = allmats_2d(idx(1:300),:);

    % out of sample components 
    outofsamp_brain = allmats_2d(idx(301:600),:); 

    [coeff,scores,~,~,explained,mu] = pca(insamp_brain);
    varexplained(:,i) = cumsum(explained);

    % multiple weights to oos brain data -> same coordinate space
    pc_brain_rep = (outofsamp_brain - mu)*coeff; 


    % insample behavior 
    tlbx_insamp_behavior = allfactors(idx(1:300),[7:13]);
    cbcl_insamp_behavior = allfactors(idx(1:300),[20:27]);
    both_insamp_behavior = allfactors(idx(1:300),[7:13 20:27]);

    % out of sample behavior 
    tlbx_outofsamp_behavior = allfactors(idx(301:600),[7:13]);
    cbcl_outofsamp_behavior = allfactors(idx(301:600),[20:27]);
    both_outofsamp_behavior = allfactors(idx(301:600),[7:13 20:27]);
    thistlbx = nan(size(scores,2),1);
    thiscbcl = nan(size(scores,2),1);
    thisboth = nan(size(scores,2),1);
    for n = 1:size(scores,2)
        % toolbox
        [A,B,~] = canoncorr(scores(:,1:n),tlbx_insamp_behavior);
        thistlbx(n,1) = corr((pc_brain_rep(:,1:n)*A(:,1)),(tlbx_outofsamp_behavior*B(:,1)));
        % toolbox
        [A,B,~] = canoncorr(scores(:,1:n),cbcl_insamp_behavior);
        thiscbcl(n,1) = corr((pc_brain_rep(:,1:n)*A(:,1)),(cbcl_outofsamp_behavior*B(:,1)));
        % toolbox
        [A,B,~] = canoncorr(scores(:,1:n),both_insamp_behavior);
        thisboth(n,1) = corr((pc_brain_rep(:,1:n)*A(:,1)),(both_outofsamp_behavior*B(:,1)));
    end
tlbx_compnum(:,i) = thistlbx;
cbcl_compnum(:,i) = thiscbcl;
both_compnum(:,i) = thisboth;
end
delete(z)

mtlbx = mean(tlbx_compnum,2);
mcbcl = mean(cbcl_compnum,2);
mboth = mean(both_compnum,2);
mve = mean(varexplained,2);
figure; hold on
h = plot(mtlbx);
h.Color = [116 198 118]./255;
h.LineWidth  = 2.5;
h = plot(mcbcl);
h.Color = [158 154 200]./255;
h.LineWidth  = 2.5;
% both
h = plot(mboth);
h.Color = [125 125 125]./255;
h.LineWidth  = 2.5;

xlim([0 300])
ylabel('Out of sample correlation')
xlabel('Number of Components')
box('off')
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/FigSx_cca_ooscorr_numcomps_n300_rsfc','tiffn')

% plot variance explained 
figure; 
plot(varexplained,'Color',[.2 .2 .2],'LineWidth',2);
xlabel('Number of Components')
ylabel('Cumulative Variance Explained (%)')
ylim([0 101])
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
box('off')
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/FigSx_cca_cumvar_numcomps_n300_rsfc','tiffn')

% save output 
cca_300_replication_numcomps.tlbxxompnum = tlbx_compnum;
cca_300_replication_numcomps.cbclcompnum = cbcl_compnum;
cca_300_replication_numcomps.bothcompnum = both_compnum;
cca_300_replication_numcomps.varexplained = varexplained;
save('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/cca_300_replication_numcomps.mat','cca_300_replication_numcomps')

%% 1000 subjects 
% initialize 1,000
tlbx_compnum = zeros(999,100);
cbcl_compnum = zeros(999,100);
both_compnum = zeros(999,100);
varexplained = zeros(999,100);
z = parpool(18);
parfor i = 1:100
    disp(['** On iteration : ' num2str(i) ' ** '])
    
    % bootstrapped index
    idx = randperm(3928);
    insamp_brain = allmats_2d(idx(1:1000),:);

    % out of sample components 
    outofsamp_brain = allmats_2d(idx(1001:2000),:); 

    % pca in-sample 
    [coeff,scores,~,~,explained,mu] = pca(insamp_brain);
    varexplained(:,i) = cumsum(explained);

    % multiple weights to oos brain data -> same coordinate space
    pc_brain_rep = (outofsamp_brain - mu)*coeff; 

    % insample behavior 
    tlbx_insamp_behavior = allfactors(idx(1:1000),[7:13]);
    cbcl_insamp_behavior = allfactors(idx(1:1000),[20:27]);
    both_insamp_behavior = allfactors(idx(1:1000),[7:13 20:27]);

    % out of sample behavior 
    tlbx_outofsamp_behavior = allfactors(idx(1001:2000),[7:13]);
    cbcl_outofsamp_behavior = allfactors(idx(1001:2000),[20:27]);
    both_outofsamp_behavior = allfactors(idx(1001:2000),[7:13 20:27]);
    thistlbx = nan(size(scores,2),1);
    thiscbcl = nan(size(scores,2),1);
    thisboth = nan(size(scores,2),1);
    for n = 1:size(scores,2)
        % toolbox
        [A,B,~] = canoncorr(scores(:,1:n),tlbx_insamp_behavior);
        thistlbx(n,1) = corr((pc_brain_rep(:,1:n)*A(:,1)),(tlbx_outofsamp_behavior*B(:,1)));
        % toolbox
        [A,B,~] = canoncorr(scores(:,1:n),cbcl_insamp_behavior);
        thiscbcl(n,1) = corr((pc_brain_rep(:,1:n)*A(:,1)),(cbcl_outofsamp_behavior*B(:,1)));
        % toolbox
        [A,B,~] = canoncorr(scores(:,1:n),both_insamp_behavior);
        thisboth(n,1) = corr((pc_brain_rep(:,1:n)*A(:,1)),(both_outofsamp_behavior*B(:,1)));
    end
tlbx_compnum(:,i) = thistlbx;
cbcl_compnum(:,i) = thiscbcl;
both_compnum(:,i) = thisboth;
end
delete(z)

mtlbx = mean(tlbx_compnum,2);
mcbcl = mean(cbcl_compnum,2);
mboth = mean(both_compnum,2);
mve = mean(varexplained,2);
figure; hold on
h = plot(mtlbx);
h.Color = [116 198 118]./255;
h.LineWidth  = 2.5;
h = plot(mcbcl);
h.Color = [158 154 200]./255;
h.LineWidth  = 2.5;
% both
h = plot(mboth);
h.Color = [125 125 125]./255;
h.LineWidth  = 2.5;

xlim([0 1000])
ylabel('Out of sample correlation')
xlabel('Number of Components')
box('off')
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/FigSx_cca_ooscorr_numcomps_n1000_rsfc','tiffn')

% plot variance explained 
figure; 
plot(varexplained,'Color',[.2 .2 .2],'LineWidth',2);
xlabel('Number of Components')
ylabel('Cumulative Variance Explained (%)')
ylim([0 101])
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
box('off')
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/FigSx_cca_cumvar_numcomps_n1000_rsfc','tiffn')

% save output 
cca_1000_replication_numcomps.tlbxxompnum = tlbx_compnum;
cca_1000_replication_numcomps.cbclcompnum = cbcl_compnum;
cca_1000_replication_numcomps.bothcompnum = both_compnum;
cca_1000_replication_numcomps.varexplained = varexplained;
save('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/cca_1000_replication_numcomps.mat','cca_1000_replication_numcomps')





%% Full sample - cortical thickness
% insample behavior 
tlbx_insamp_behavior = allfactors(1:1814,7:13);
cbcl_insamp_behavior = allfactors(1:1814,[20:27]);
both_insamp_behavior = allfactors(1:1814,[7:13 20:27]);

% out of sample behavior 
tlbx_outofsamp_behavior = allfactors(1815:end,[7:13]);
cbcl_outofsamp_behavior = allfactors(1815:end,[20:27]);
both_outofsamp_behavior = allfactors(1815:end,[7:13 20:27]);

tlbx_compnum = zeros(1813,100);
cbcl_compnum = zeros(1813,100);
both_compnum = zeros(1813,100);
varexplained = zeros(1813,100);
z = parpool(18);
parfor i = 1:100
    disp(['** On iteration : ' num2str(i) ' ** '])
    
    % bootstrapped index
    idx = randperm(3604);
    insamp_brain = allmats_2d(idx(1:1814),:);

    % out of sample components 
    outofsamp_brain = allmats_2d(idx(1815:3604),:); 

    % pca in-sample 
    [coeff,scores,~,~,explained,mu] = pca(insamp_brain);
    varexplained(:,i) = cumsum(explained);

    % multiple weights to oos brain data -> same coordinate space
    pc_brain_rep = (outofsamp_brain - mu)*coeff; 


    % insample behavior 
    tlbx_insamp_behavior = allfactors(idx(1:1814),[7:13]);
    cbcl_insamp_behavior = allfactors(idx(1:1814),[20:27]);
    both_insamp_behavior = allfactors(idx(1:1814),[7:13 20:27]);

    % out of sample behavior 
    tlbx_outofsamp_behavior = allfactors(idx(1815:end),[7:13]);
    cbcl_outofsamp_behavior = allfactors(idx(1815:end),[20:27]);
    both_outofsamp_behavior = allfactors(idx(1815:end),[7:13 20:27]);
    thistlbx = nan(size(scores,2),1);
    thiscbcl = nan(size(scores,2),1);
    thisboth = nan(size(scores,2),1);
    for n = 1:size(scores,2)
        % toolbox
        [A,B,~] = canoncorr(scores(:,1:n),tlbx_insamp_behavior);
        thistlbx(n,1) = corr((pc_brain_rep(:,1:n)*A(:,1)),(tlbx_outofsamp_behavior*B(:,1)));
        % toolbox
        [A,B,~] = canoncorr(scores(:,1:n),cbcl_insamp_behavior);
        thiscbcl(n,1) = corr((pc_brain_rep(:,1:n)*A(:,1)),(cbcl_outofsamp_behavior*B(:,1)));
        % toolbox
        [A,B,~] = canoncorr(scores(:,1:n),both_insamp_behavior);
        thisboth(n,1) = corr((pc_brain_rep(:,1:n)*A(:,1)),(both_outofsamp_behavior*B(:,1)));
    end
tlbx_compnum(:,i) = thistlbx;
cbcl_compnum(:,i) = thiscbcl;
both_compnum(:,i) = thisboth;
end
delete(z)

mtlbx = mean(tlbx_compnum,2);
mcbcl = mean(cbcl_compnum,2);
mboth = mean(both_compnum,2);
mve = mean(varexplained,2);
figure; hold on
h = plot(mtlbx);
h.Color = [116 198 118]./255;
h.LineWidth  = 2.5;
h = plot(mcbcl);
h.Color = [158 154 200]./255;
h.LineWidth  = 2.5;
% both
h = plot(mboth);
h.Color = [125 125 125]./255;
h.LineWidth  = 2.5;

xlim([0 1965])
ylabel('Out of sample correlation')
xlabel('Number of Components')
box('off')
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/FigSx_cca_ooscorr_numcomps_nfs_thickness','tiffn')

% plot variance explained 
figure; 
plot(varexplained,'Color',[.2 .2 .2],'LineWidth',2);
xlabel('Number of Components')
ylabel('Cumulative Variance Explained (%)')
ylim([0 101])
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/FigSx_cca_cumvar_numcomps_nfs_thickness','tiffn')

% save output 
cca_fs_replication_numcomps.tlbxcompnum = tlbx_compnum;
cca_fs_replication_numcomps.cbclcompnum = cbcl_compnum;
cca_fs_replication_numcomps.bothcompnum = both_compnum;
cca_fs_replication_numcomps.varexplained = varexplained;
save('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/cca_fs_replication_numcomps_thickness.mat','cca_fs_replication_numcomps')

%% 300 subjects 
% initialize 300
tlbx_compnum = zeros(299,100);
cbcl_compnum = zeros(299,100);
both_compnum = zeros(299,100);
varexplained = zeros(299,100);
z = parpool(18);
parfor i = 1:100
    disp(['** On iteration : ' num2str(i) ' ** '])
    
    % bootstrapped index
    idx = randperm(3604);
    insamp_brain = allmats_2d(idx(1:300),:);

    % out of sample components 
    outofsamp_brain = allmats_2d(idx(301:600),:); 

    % pca in-sample 
    [coeff,scores,~,~,explained,mu] = pca(insamp_brain);
    varexplained(:,i) = cumsum(explained);

    % multiple weights to oos brain data -> same coordinate space
    pc_brain_rep = (outofsamp_brain - mu)*coeff; 

    % insample behavior 
    tlbx_insamp_behavior = allfactors(idx(1:300),[7:13]);
    cbcl_insamp_behavior = allfactors(idx(1:300),[20:27]);
    both_insamp_behavior = allfactors(idx(1:300),[7:13 20:27]);

    % out of sample behavior 
    tlbx_outofsamp_behavior = allfactors(idx(301:600),[7:13]);
    cbcl_outofsamp_behavior = allfactors(idx(301:600),[20:27]);
    both_outofsamp_behavior = allfactors(idx(301:600),[7:13 20:27]);
    thistlbx = nan(size(scores,2),1);
    thiscbcl = nan(size(scores,2),1);
    thisboth = nan(size(scores,2),1);
    for n = 1:size(scores,2)
        % toolbox
        [A,B,~] = canoncorr(scores(:,1:n),tlbx_insamp_behavior);
        thistlbx(n,1) = corr((pc_brain_rep(:,1:n)*A(:,1)),(tlbx_outofsamp_behavior*B(:,1)));
        % toolbox
        [A,B,~] = canoncorr(scores(:,1:n),cbcl_insamp_behavior);
        thiscbcl(n,1) = corr((pc_brain_rep(:,1:n)*A(:,1)),(cbcl_outofsamp_behavior*B(:,1)));
        % toolbox
        [A,B,~] = canoncorr(scores(:,1:n),both_insamp_behavior);
        thisboth(n,1) = corr((pc_brain_rep(:,1:n)*A(:,1)),(both_outofsamp_behavior*B(:,1)));
    end
tlbx_compnum(:,i) = thistlbx;
cbcl_compnum(:,i) = thiscbcl;
both_compnum(:,i) = thisboth;
end
delete(z)

mtlbx = mean(tlbx_compnum,2);
mcbcl = mean(cbcl_compnum,2);
mboth = mean(both_compnum,2);
mve = mean(varexplained,2);
figure; hold on
h = plot(mtlbx);
h.Color = [116 198 118]./255;
h.LineWidth  = 2.5;
h = plot(mcbcl);
h.Color = [158 154 200]./255;
h.LineWidth  = 2.5;
% both
h = plot(mboth);
h.Color = [125 125 125]./255;
h.LineWidth  = 2.5;

xlim([0 300])
ylabel('Out of sample correlation')
xlabel('Number of Components')
box('off')
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/FigSx_cca_ooscorr_numcomps_n300_thickness','tiffn')

% plot variance explained 
figure; 
plot(varexplained,'Color',[.2 .2 .2],'LineWidth',2);
xlabel('Number of Components')
ylabel('Cumulative Variance Explained (%)')
ylim([0 101])
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/FigSx_cca_cumvar_numcomps_n300_thickness','tiffn')

% save output 
cca_300_replication_numcomps.tlbxxompnum = tlbx_compnum;
cca_300_replication_numcomps.cbclcompnum = cbcl_compnum;
cca_300_replication_numcomps.bothcompnum = both_compnum;
cca_300_replication_numcomps.varexplained = varexplained;
save('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/cca_300_replication_numcomps_thickness.mat','cca_300_replication_numcomps')

%% 1000 subjects 
% initialize 1,000
tlbx_compnum = zeros(999,100);
cbcl_compnum = zeros(999,100);
both_compnum = zeros(999,100);
varexplained = zeros(999,100);
z = parpool(18);
parfor i = 1:100
    disp(['** On iteration : ' num2str(i) ' ** '])
    
    % bootstrapped index
    idx = randperm(3604);
    insamp_brain = allmats_2d(idx(1:1000),:);

    % out of sample components 
    outofsamp_brain = allmats_2d(idx(1001:2000),:); 

    % pca in-sample 
    [coeff,scores,~,~,explained,mu] = pca(insamp_brain);
    varexplained(:,i) = cumsum(explained);

    % multiple weights to oos brain data -> same coordinate space
    pc_brain_rep = (outofsamp_brain - mu)*coeff; 


    % insample behavior 
    tlbx_insamp_behavior = allfactors(idx(1:1000),[7:13]);
    cbcl_insamp_behavior = allfactors(idx(1:1000),[20:27]);
    both_insamp_behavior = allfactors(idx(1:1000),[7:13 20:27]);

    % out of sample behavior 
    tlbx_outofsamp_behavior = allfactors(idx(1001:2000),[7:13]);
    cbcl_outofsamp_behavior = allfactors(idx(1001:2000),[20:27]);
    both_outofsamp_behavior = allfactors(idx(1001:2000),[7:13 20:27]);
    thistlbx = nan(size(scores,2),1);
    thiscbcl = nan(size(scores,2),1);
    thisboth = nan(size(scores,2),1);
    for n = 1:size(scores,2)
        % toolbox
        [A,B,~] = canoncorr(scores(:,1:n),tlbx_insamp_behavior);
        thistlbx(n,1) = corr((pc_brain_rep(:,1:n)*A(:,1)),(tlbx_outofsamp_behavior*B(:,1)));
        % toolbox
        [A,B,~] = canoncorr(scores(:,1:n),cbcl_insamp_behavior);
        thiscbcl(n,1) = corr((pc_brain_rep(:,1:n)*A(:,1)),(cbcl_outofsamp_behavior*B(:,1)));
        % toolbox
        [A,B,~] = canoncorr(scores(:,1:n),both_insamp_behavior);
        thisboth(n,1) = corr((pc_brain_rep(:,1:n)*A(:,1)),(both_outofsamp_behavior*B(:,1)));
    end
tlbx_compnum(:,i) = thistlbx;
cbcl_compnum(:,i) = thiscbcl;
both_compnum(:,i) = thisboth;
end
delete(z)

mtlbx = mean(tlbx_compnum,2);
mcbcl = mean(cbcl_compnum,2);
mboth = mean(both_compnum,2);
mve = mean(varexplained,2);
figure; hold on
h = plot(mtlbx);
h.Color = [116 198 118]./255;
h.LineWidth  = 2.5;
h = plot(mcbcl);
h.Color = [158 154 200]./255;
h.LineWidth  = 2.5;
% both
h = plot(mboth);
h.Color = [125 125 125]./255;
h.LineWidth  = 2.5;

xlim([0 1000])
ylabel('Out of sample correlation')
xlabel('Number of Components')
box('off')
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/FigSx_cca_ooscorr_numcomps_n1000_thickness','tiffn')

% plot variance explained 
figure; 
plot(varexplained,'Color',[.2 .2 .2],'LineWidth',2);
xlabel('Number of Components')
ylabel('Cumulative Variance Explained (%)')
ylim([0 101])
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/FigSx_cca_cumvar_numcomps_n1000_thickness','tiffn')

% save output 
cca_1000_replication_numcomps.tlbxxompnum = tlbx_compnum;
cca_1000_replication_numcomps.cbclcompnum = cbcl_compnum;
cca_1000_replication_numcomps.bothcompnum = both_compnum;
cca_1000_replication_numcomps.varexplained = varexplained;
save('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/cca_1000_replication_numcomps_thickness.mat','cca_1000_replication_numcomps')
