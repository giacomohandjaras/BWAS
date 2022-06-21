% make reliability histograms for UKB, ABCD, HCP

hcp = load('/data/nil-bluearc/GMT/Scott/HCP/Parcelcorrmats/Splithalf_reliability.mat');
hcp = hcp.Reliability;

abcd_reliability = load('/data/nil-bluearc/GMT/Scott/ABCD/Parcelcorrmats/SplitHalfReliability_rest.mat');
abcd_reliability_rep = load('/data/nil-bluearc/GMT/Scott/ABCD/Parcelcorrmats/SplitHalfReliability_rest_rep.mat');
abcd = [abcd_reliability.Reliability ; abcd_reliability_rep.Reliability];

ukb = load('/data/nil-bluearc/GMT/Scott/UKB/rsfc_reliability.mat');
ukb = ukb.ukb_reliability;

figure; hold on 

% UKB
h=histogram(ukb);
h.FaceColor = [241 127 116]./255;

i = histogram(abcd);
i.FaceColor = [175 97 123]./255;

j = histogram(hcp(hcp>0));
j.FaceColor = [98 76 105]./255;

saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/ukb_abcd_hcp_rsfc_reliability','tiffn')
