%% Figure 4 

% Probability (0 - 1) of Type 1 error, type 2 error, inflation, sign,
% multivariate replication 

abcd_make_big_data_table_FD_08

abcd_make_brainbehavior_variables

%% Sampling variability 
load('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/EffectSizeSimilarity.mat');

for f = 1:size(EffectSizeSimilarity,3)
   this = EffectSizeSimilarity(:,:,f); 
   mm(:,f) = nanmean(this,2);
end
mm = nanmean(this,2);  
    
    
%mthisbin = mean(Thesecorrs,1);
%mm = sort((max(Thesecorrs) - min(Thesecorrs)) ./ 2,'descend');
%mm = (mm-min(mm))/(max(mm)-min(mm));

%% Type 2 error 

load('/data/nil-bluearc/GMT/Scott/ABCD/type2error_41vars_rsfc_pval.mat');
t2 = squeeze(nanmean(type2error,3));%./100;
t2 = t2(:,end)./100;
% Normalize 
%t2 = (t2-min(t2))/(max(t2)-min(t2));

%% Type 1 error 

%load('/data/nil-bluearc/GMT/Scott/ABCD/type1error_41vars_rsfc.mat');
%t1 = squeeze(type1error(:,4,17))./100;
% Normalize 
%t1 = (t1-min(t1))/(max(t1)-min(t1));

%% Inflation 
load('/data/nil-bluearc/GMT/Scott/ABCD/typemerror_41vars_rsfc_pvals.mat');
mic = [];
for f = 1:size(typemerror_pval,1)
    this = typemerror_pval{f};
    this = this{7};
    % 50% inflation 
    mic(:,f) = mean(squeeze(this(:,:,4)),1);
end
mic = nanmean(mic,2);
mic(isnan(mic)) = 100;
mic=mic./100;
%mic = (mic-min(mic))/(max(mic)-min(mic));


%% Sign flip 
load('/data/nil-bluearc/GMT/Scott/ABCD/typeserror_41vars_rsfc_pval.mat');

msf = [];
for f = 1:size(typeserror_pval,1)
    this = typeserror_pval{f};
    this = this{7};
    % 50% inflation 
    msf(:,f) = mean(this,1);
end
msf = nanmean(msf,2);
msf = msf./100;

%% sVR 
% RSFC
%IQ 
rsfc = readtable('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/PCA5update20200819.rsfc.csv' , 'Delimiter',',');
rsfc.traincor(3200) = .5;
rsfc.outofsamplecor(3200) = .34;
discovery = reshape(rsfc.traincor(1601:3200),100,16)';
replication = reshape(rsfc.outofsamplecor(1601:3200),100,16)';

thismean = mean(discovery,2);
maximumd = thismean(end);

meanr = mean(replication,2) ./ maximumd;

% both 
%cca = mean(outofsampleweightsall.both,2);
%cca = sort((cca-min(cca))/(max(cca)-min(cca)));

%% Plot everything 
cumprob = t2 + mic + msf;
cumprob(cumprob>1) = 1;
%cumprob = (cumprob-min(cumprob))/(max(cumprob)-min(cumprob));

figure; hold on 
% cumulative probability of comitting a statistical error (1,2,m,s)
plot(1:16,cumprob,'Color',[25 0 0]./255,'LineWidth',1);
% cca 
plot(1:16,meanr,'Color',[0 100 0]./255,'LineWidth',1);
% sv 
plot(1:16,mm.*-1,'Color',[0 0 150]./255,'LineWidth',1);

xlim([1 16])
ylim([0 1])
xlabel('Number of Subjects','FontWeight','bold','FontSize',14)
ylabel('Normalized Replication Coefficient','FontWeight','bold','FontSize',14)
title('RSFC')
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/Normalized_all_fig5','tiffn')





