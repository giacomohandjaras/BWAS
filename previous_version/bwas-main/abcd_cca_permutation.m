%% CCA permuted null: full sample, RSFC & Cortical thickness 

abcd_make_big_data_table_FD_08

abcd_make_brainbehavior_variables


%% Run Model - RSFC 

% Set input params 
iter = 100; % null interations

% brain data 
insamp_brain = allmats_2d(1:1964,:); 
outofsamp_brain = allmats_2d(1965:end,:);

% PCA on in sample; apply coordinates to out of sample 
[coeffs,scores,~,~,explained,mu] = pca(insamp_brain);
pc_brain_rep = (outofsamp_brain - mu)*coeffs; % Put left out set in same coordinate space 

% Get 20% of variance explained 
numcomponents = find(cumsum(explained)>=20);
numcomponents = numcomponents(1);

% behavioral data - toolbox 
insamp_behavior = allfactors(1:size(insamp_brain,1),[7:13]);
outofsamp_behavior = allfactors(1+size(insamp_brain,1):end,[7:13]);

% Permutation
theseweights = zeros(1,iter);
inweights = zeros(1,iter);
for j = 1:iter  
    idx = randperm(size(insamp_brain,1));
     % Grab index where cumulative var explained reaches 20% 
    [A,B,r] = canoncorr(scores(:,1:numcomponents),insamp_behavior(idx,:));
    inweights(1,j) = r(1);
    theseweights(1,j) = corr((pc_brain_rep(:,1:numcomponents)*A(:,1)),(outofsamp_behavior*B(:,1)));
end
% Output 
rsfc_insampleweights_tlbx = inweights;
rsfc_outofsampleweights_tlbx = theseweights;

% behavioral data - CBCL
insamp_behavior = allfactors(1:size(insamp_brain,1),[20:27]);
outofsamp_behavior = allfactors(1+size(insamp_brain,1):end,[20:27]);

% Permutation
theseweights = zeros(1,iter);
inweights = zeros(1,iter);
for j = 1:iter  
    idx = randperm(size(insamp_brain,1));
     % Grab index where cumulative var explained reaches 50% 
    [A,B,r] = canoncorr(scores(:,1:numcomponents),insamp_behavior(idx,:));
    inweights(1,j) = r(1);
    theseweights(1,j) = corr((pc_brain_rep(:,1:numcomponents)*A(:,1)),(outofsamp_behavior*B(:,1)));
end
% Output 
rsfc_insampleweights_cbcl = inweights;
rsfc_outofsampleweights_cbcl = theseweights;

% % behavioral data - CBCL & Toolbox 
% insamp_behavior = allfactors(1:size(insamp_brain,2),[2:8 12:18]);
% outofsamp_behavior = allfactors(1+size(insamp_brain,2):end,[2:8 12:18]);
% 
% % Permutation
% theseweights = zeros(1,iter);
% inweights = zeros(1,iter);
% for j = 1:iter  
%     idx = randperm(size(insamp_brain,2));
%      % Grab index where cumulative var explained reaches 50% 
%     [A,B,r] = canoncorr(V(:,1:numcomponents),insamp_behavior(idx,:));
%     inweights(1,j) = r(1);
%     theseweights(1,j) = corr((ooscorrds(:,1:numcomponents)*A(:,1)),(outofsamp_behavior*B(:,1)));
% end
% % Output 
% rsfc_insampleweights_both = inweights;
% rsfc_outofsampleweights_both = theseweights;

%rsfc_permuted_cca.rsfc_insampleweights_both = rsfc_insampleweights_both;
rsfc_permuted_cca.rsfc_insampleweights_tlbx = rsfc_insampleweights_tlbx;
rsfc_permuted_cca.rsfc_insampleweights_cbcl = rsfc_insampleweights_cbcl;
%rsfc_permuted_cca.rsfc_outofsampleweights_both = rsfc_outofsampleweights_both;
rsfc_permuted_cca.rsfc_outofsampleweights_tlbx = rsfc_outofsampleweights_tlbx;
rsfc_permuted_cca.rsfc_outofsampleweights_cbcl = rsfc_outofsampleweights_cbcl;

% Save
save('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/CCA/rsfc_permuted_cca.mat','rsfc_permuted_cca');

%% Run model - cortical thickness 
clear all 
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
iter = 1000;

% brain
insamp_brain = allmats_2d(1:1814,:); 
outofsamp_brain = allmats_2d(1815:end,:); 

[coeffs,scores,~,~,explained,mu] = pca(insamp_brain);
% PCA on in sample; apply coordinates to out of sample 
pc_brain_rep = (outofsamp_brain - mu)*coeffs; % Put left out set in same coordinate space

% Get 20% of variance explained 
numcomponents = find(cumsum(explained)>=20);
numcomponents = numcomponents(1);

% Toolbox
insamp_behavior = allfactors(1:1814,[7:13]);
outofsamp_behavior = allfactors(1815:end,[7:13]);

% Permutation
theseweights = zeros(1,iter);
inweights = zeros(1,iter);
for j = 1:iter  
    idx = randperm(size(insamp_brain,1));
     % Grab index where cumulative var explained reaches 50% 
    [A,B,r] = canoncorr(scores(:,1:numcomponents),insamp_behavior(idx,:));
    inweights(1,j) = r(1);
    theseweights(1,j) = corr((pc_brain_rep(:,1:numcomponents)*A(:,1)),(outofsamp_behavior*B(:,1)));
end
% Output 
thickness_insampleweights_tlbx = inweights;
thickness_outofsampleweights_tlbx = theseweights;

% behavioral data - CBCL
insamp_behavior = allfactors(1:size(insamp_brain,1),[20:27]);
outofsamp_behavior = allfactors(1+size(insamp_brain,1):end,[20:27]);

% Permutation
theseweights = zeros(1,iter);
inweights = zeros(1,iter);
for j = 1:iter  
    idx = randperm(size(insamp_brain,1));
     % Grab index where cumulative var explained reaches 50% 
    [A,B,r] = canoncorr(scores(:,1:numcomponents),insamp_behavior(idx,:));
    inweights(1,j) = r(1);
    theseweights(1,j) = corr((pc_brain_rep(:,1:numcomponents)*A(:,1)),(outofsamp_behavior*B(:,1)));
end
% Output 
thickness_insampleweights_cbcl = inweights;
thickness_outofsampleweights_cbcl = theseweights;


thickness_permuted_cca.thickness_insampleweights_tlbx = thickness_insampleweights_tlbx;
thickness_permuted_cca.thickness_insampleweights_cbcl = thickness_insampleweights_cbcl;
thickness_permuted_cca.thickness_outofsampleweights_tlbx = thickness_outofsampleweights_tlbx;
thickness_permuted_cca.thickness_outofsampleweights_cbcl = thickness_outofsampleweights_cbcl;

% Save
save('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/CCA/thickness_permuted_cca.mat','thickness_permuted_cca');


%% Get 95th percentile 

% RSFC
load('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/CCA/rsfc_permuted_cca.mat')

% percentiles

% toolbox 
ninetyninthprctile_tlbx_rsfc = prctile(rsfc_permuted_cca.rsfc_outofsampleweights_tlbx,99);

% cbcl
ninetyninthprctile_cbcl_rsfc = prctile(rsfc_permuted_cca.rsfc_outofsampleweights_cbcl,99);

% both 
ninetyninthprctile_both_rsfc = prctile(rsfc_permuted_cca.rsfc_outofsampleweights_both,99);


% Cortical Thickness
load('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/CCA/thickness_permuted_cca.mat')

% percentiles

% toolbox 
ninetyninthprctile_tlbx_thickness = prctile(thickness_permuted_cca.thickness_outofsampleweights_tlbx,99);

% cbcl
ninetyninthprctile_cbcl_thickness = prctile(thickness_permuted_cca.thickness_outofsampleweights_cbcl,99);

% both 
ninetyninthprctile_both_thickness = prctile(thickness_permuted_cca.thickness_outofsampleweights_both,99);


% Compare to observed values 
% RSFC
load('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/CCA/outofsampleweightsall_50percentcomps.mat')

% Toolbox 
rsfc_observed = mean(outofsampleweightsall.tlbx,2);
rsfc = rsfc_observed(end); % 1st column is observed 
rsfc(1,2) =  ninetyfifthprctile_tlbx_rsfc; % 2nd column is permuted 99th percentile 

% CBCL 
rsfc_observed = mean(outofsampleweightsall.cbcl,2);
rsfc(2,1) = rsfc_observed(end); % 1st column is observed 
rsfc(2,2) =  ninetyfifthprctile_cbcl_rsfc; % 2nd column is permuted 99th percentile 

% Both
rsfc_observed = mean(outofsampleweightsall.both,2);
rsfc(3,1) = rsfc_observed(end); % 1st column is observed 
rsfc(3,2) =  ninetyfifthprctile_both_rsfc; % 2nd column is permuted 99th percentile 


% thickness
load('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/CCA/outofsampleweightsall_50percentcomp_thickness.mat')

% Toolbox 
ct_observed = mean(outofsampleweightsall.tlbx,2);
ct = ct_observed(end); % 1st column is observed 
ct(1,2) =  ninetyfifthprctile_tlbx_thickness; % 2nd column is permuted 99th percentile 

% CBCL 
ct_observed = mean(outofsampleweightsall.cbcl,2);
ct(2,1) = ct_observed(end); % 1st column is observed 
ct(2,2) =  ninetyfifthprctile_cbcl_thickness; % 2nd column is permuted 99th percentile 

% Both
ct_observed = mean(outofsampleweightsall.both,2);
ct(3,1) = ct_observed(end); % 1st column is observed 
ct(3,2) =  ninetyfifthprctile_both_thickness; % 2nd column is permuted 99th percentile 

disp('RSFC values: ')
rsfc
disp('CT values: ')
ct

