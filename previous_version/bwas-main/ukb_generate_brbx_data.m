addpath(genpath('/data/cn/data1/scripts/CIFTI_RELATED/'))

datadir = '/data/Daenerys/UKB/ukbdownload/';

oldbehavioral = readtable([datadir 'behavioral_metrics_20210831.csv']);
behavioral = readtable([datadir 'behavioral_data']);

key = readtable([datadir 'behavioral_metrics_20210831_key.txt']);

subjecttable = innerjoin(behavioral,oldbehavioral,'Keys','eid');
subjecttable(:,3:end) = [];
id = {'NaN' '' 'NA' 999 777 NaN Inf '999' '777'};
TF = ismissing(subjecttable,id);

% Remove subjects with missing data 
subjecttable(TF(:,2)>0,:) = [];

for s = 1:size(subjecttable,1)
    FIQ(s,1) = str2num(subjecttable.x20016_2_0{s});
end
subjecttable = addvars(subjecttable,FIQ);
subjecttable(:,2) = [];
subjects = subjecttable.eid;

%% Loop thru pconns 
allmats2d = nan(55278,size(subjects,1));
Missingidx = zeros(size(subjects,1),1);
z=parpool(15);
parfor s = 1:size(subjects,1)
    disp(['Subject ' num2str(s)])
    if exist([datadir 'sub-' num2str(subjects(s,1)) '/ses-01/sub-' num2str(subjects(s,1)) '_ses-01_task-rest_bold_atlas-Gordon2014FreeSurferSubcortical_desc-filtered_timeseries.ptseries.nii_all_frames_at_FD_0.2.pconn.nii'],'file')
        thispconn = ft_read_cifti_mod([datadir 'sub-' num2str(subjects(s,1)) '/ses-01/sub-' num2str(subjects(s,1)) '_ses-01_task-rest_bold_atlas-Gordon2014FreeSurferSubcortical_desc-filtered_timeseries.ptseries.nii_all_frames_at_FD_0.2.pconn.nii']);
        rmat = thispconn.data(1:333,1:333);
        uidx = find(triu(rmat,1));
        allmats2d(:,s) = rmat(uidx);
    else
        disp(['Missing subject ' num2str(s)])
        Missingidx(s,1) = 1;
    end
    
end
delete(z)

%% Behavioral data formatting 

% subjecttable = table(subjects);
% behavioral.Properties.VariableNames{2} = 'subjects';
% behavioraltable = innerjoin(behavioral,subjecttable,'Keys','subjects');
% behavioraltable = removevars(behavioraltable, 'x20157_0_0');
% behavioraltable = removevars(behavioraltable, {'x20240_0_0','x20433_0_0','x20434_0_0','x20442_0_0','x20458_0_0','x20506_0_0'});
% 
% id = {'NaN' '' 999 777 NaN Inf '999' '777'};
% TF = ismissing(behavioraltable,id);
% 
% TFs = sum(TF')';
% 
% % Remove subjects with missing behavioral data 
% behavioraltable(TFs>0,:) = [];
mats = allmats2d'; clear allmats2d

% Remove subjects with missing rsfc data
a=find(sum(isnan(mats'))>0)';
mats(a,:) = [];
subjecttable(a,:) = [];

for i = 1:size(mats,2)
    ukbcorrs(i,1) = corr(mats(:,i),subjecttable.FIQ);
end

%% subsample to N=877, 100 times 
ukbsubsample = nan(size(mats,2),100);
z = parpool(15);
parfor i = 1:100
    idx = datasample(1:size(mats,1),900)';
    this = [];
    for r = 1:size(mats,2)
       this(r,1) = corr(mats(idx,r),subjecttable.FIQ(idx)); 
    end
    ukbsubsample(:,i) = this;
end
delete(z)

%% HCP
% Load HCP 
HCP_make_big_datatable
%
rmats = cell2mat(MainTable_HCP.Corrmats);
rmats = reshape(rmats,333,size(MainTable_HCP,1),333);
rmats = permute(rmats,[1,3,2]);

% Behavioral measure(s) (subject x measure)
Behavioralitems = MainTable_HCP.CogFluidComp_Unadj(MainTable_HCP.FrameTotal >= 667 & MainTable_HCP.MissingRestSubjects == 0 & MainTable_HCP.Badtmaskidx == 0);
missing = isnan(Behavioralitems);
nanidx = find(missing==0);
Behavioralitems = Behavioralitems(nanidx);

% All correlation matrices, filtered and ordered
rmats = rmats(:,:,MainTable_HCP.FrameTotal >= 667 & MainTable_HCP.MissingRestSubjects == 0 & MainTable_HCP.Badtmaskidx == 0);
rmats = rmats(:,:,nanidx);

% All edges (subject x edges)
uidx = find(triu(rmats(:,:,1),1));
rmats_2d = zeros(size(rmats,3),length(uidx));
for s = 1:size(rmats,3)
    Thismat = rmats(:,:,s);
    rmats_2d(s,:) = Thismat(uidx);
end

% Correlations with behavior 
hcpsubsample = zeros(size(rmats_2d,2),size(Behavioralitems,2));
for i = 1:size(Behavioralitems,2)
    disp(['On variable ' num2str(i)])
    ThisBehavior = Behavioralitems(:,i);
    for x = 1:size(hcpsubsample,1)
        hcpsubsample(x,i) = corr(rmats_2d(:,x),ThisBehavior);
    end
end

%% ABCD 
abcd_make_big_data_table_FD_08
abcd_make_brainbehavior_variables

abcdmats = allmats(1:333,1:333,:);
allmats_2d = [];

% All edges (subject x edges)
uidx = find(triu(abcdmats(:,:,1),1));
allmats_2d = zeros(size(abcdmats,3),length(uidx));
for s = 1:size(abcdmats,3)
    Thismat = abcdmats(:,:,s);
    allmats_2d(s,:) = Thismat(uidx);
end

vars = allfactors(:,14);

% Full sample 
abcdcorrs = [];
for r = 1:size(allmats_2d,2)
   abcdcorrs(r,1) = corr(allmats_2d(:,r),vars); 
end

% subsampled
abcdr = zeros(size(allmats_2d,2),100);
z = parpool(15);
parfor i = 1:100
    idx = datasample(1:size(allmats_2d,1),900);
    thisvar = [];
    for r = 1:size(allmats_2d,2)
        thisvar(r,1) = corr(allmats_2d(idx,r),vars(idx));
    end
    abcdr(:,i) = thisvar;    
end
abcdsubsample = abcdr;
delete(z)

% Compare percentiles in fluid IQ 
disp(['Uk Biobank: ' num2str(prctile(abs(ukbsubsample(:)),99))])
disp(['ABCD: ' num2str(prctile(abs(abcdsubsample(:)),99))])
disp(['HCP: ' num2str(prctile(abs(hcpsubsample),99))])


%% Make histograms 

% Subsample
figure;
% Full sample  (B)
subplot(1,2,2);
hold on

% HCP
j = histfit(hcpsubsample);
j(1).YData = [];
j(2).Color = [98 76 105]./255;

% ABCD
i = histfit(abcdcorrs);
i(1).YData = [];
i(2).Color = [175 97 123]./255;

% UKB
h=histfit(ukbcorrs);
h(1).YData = [];
h(2).Color = [241 127 116]./255;

xlim([-.2 .2])
xlabel('Correlation (r)')
ylabel('Count')
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)

subplot(1,2,1); 
hold on 


%hcp
h1=histfit(hcpsubsample);
h1(1).YData = [];
h1(2).Color = [98 76 105]./255;
% h1.Normalization = 'probability';
% h1.BinWidth = 0.01;
% h1.FaceColor = [98 76 105]./255;
% h1.EdgeColor = [98 76 105]./255;
% h1.FaceAlpha = fa;
% h1.EdgeAlpha = ea;

% Get 55k random samples and plot 
idx=randperm(length(abcdsubsample(:)))';
abcds = abcdsubsample(:);
h2 = histfit(abcds(idx(1:size(hcpsubsample,1))));
h2(1).YData = [];
h2(2).Color = [175 97 123]./255;
% h2.Normalization = 'probability';
% h2.BinWidth = 0.01;
% h2.FaceColor = [107 174 214]./255;
% h2.EdgeColor = [107 174 214]./255;
% h2.FaceAlpha = fa;
% h2.EdgeAlpha = ea;

idx=randperm(length(ukbsubsample(:)))';
ukbs = ukbsubsample(:);
h3 = histfit(ukbs(idx(1:size(hcpsubsample,1))));
h3(1).YData = [];
h3(2).Color = [241 127 116]./255;
% h3 = histfit(ukbsubsample(:));
% h3.Normalization = 'probability';
% h3.BinWidth = 0.01;
% h3.FaceColor = [198 219 239]./255;
% h3.EdgeColor = [198 219 239]./255;
% h3.FaceAlpha = fa;
% h3.EdgeAlpha = ea;

xlim([-.2 .2])
ylim([0 1000])
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
box('off')


