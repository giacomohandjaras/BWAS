%% plot replication weights as a function of discovery weights, across sample sizes
% for RSFC (tlbx, cbcl, nih-t-cbcl) and thickness (tlbx, cbcl, nih-t-cbcl)

%% RSFC 

% load in sample and out of sample weight s
load('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/CCA/insampleweightsall_20percentcomps.mat');
load('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/CCA/outofsampleweightsall_20percentcomps.mat');

% plot Toolbox 
discovery = insampleweightsall.tlbx;
replication = outofsampleweightsall.tlbx;

% color pallete
cmap = [];
cmap(1,:) = [0 68 27]./255;
cmap(2,:) = [65 171 93]./255;
cmap(3,:) = [199 233 192]./255;
[x,y]=meshgrid([1:size(cmap,1)],1:16);
cmap = interp2(x([1,8,16],:),y([1,8,16],:),cmap,x,y);

% plot 
figure; 
subplot(2,4,7)
hold on 
for s = 1:size(insampleweightsall.tlbx,1)
    h = scatter(replication(s,:),discovery(s,:));
    h.Marker = '.';
    h.SizeData = 50; 
    h.CData = cmap(s,:);
end
xlim([-.15 .4])
ylim([.2 1])
h = scatter(mean(replication(16,:)),mean(discovery(16,:)));
h.Marker = '.';
h.SizeData = 1000; 
h.CData = [0 0 0];
h = line([.094 .094],[.2 1]);
h.LineStyle = '--';
h.LineWidth = 2;
h.Color = [.4 .4 .4];
l = scatter(mean(replication(8,:)),mean(discovery(8,:)));
l.Marker = '.';
l.SizeData = 1000; 
l.CData = [.6 0 0];

% colorbar 
%figure; scatter(discovery(:,1),replication(:,1),20,cmap,'filled');colormap(cmap);colorbar
%saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/Fig2A_C_colorbar','tiffn')

% plot CBCL 
discovery = insampleweightsall.cbcl;
replication = outofsampleweightsall.cbcl;

% color pallete
cmap = [];
cmap(1,:) = [63 0 125]./255;
cmap(2,:) = [128 125 186]./255;
cmap(3,:) = [218 218 235]./255;
[x,y]=meshgrid([1:size(cmap,1)],1:16);
cmap = interp2(x([1,8,16],:),y([1,8,16],:),cmap,x,y);

% plot 
subplot(2,4,8)
hold on 
for s = 1:size(insampleweightsall.cbcl,1)
    h = scatter(replication(s,:),discovery(s,:));
    h.Marker = '.';
    h.SizeData = 50; 
    h.CData = cmap(s,:);
end
xlim([-.15 .4])
ylim([.2 1])
h = scatter(mean(replication(16,:)),mean(discovery(16,:)));
h.Marker = '.';
h.SizeData = 1000; 
h.CData = [0 0 0];
h = line([.07 .07],[.2 1]);
h.LineStyle = '--';
h.LineWidth = 2;
h.Color = [.4 .4 .4];
l = scatter(mean(replication(8,:)),mean(discovery(8,:)));
l.Marker = '.';
l.SizeData = 1000; 
l.CData = [.6 0 0];


% % plot NIH-T-CBCL 
% discovery = insampleweightsall.both;
% replication = outofsampleweightsall.both;
% 
% % color pallete
% cmap = [];
% cmap(1,:) = [0 0 0]./255;
% cmap(2,:) = [115 115 115]./255;
% cmap(3,:) = [217 217 217]./255;
% [x,y]=meshgrid([1:size(cmap,1)],1:16);
% cmap = interp2(x([1,8,16],:),y([1,8,16],:),cmap,x,y);
% 
% % plot 
% subplot(1,3,3)
% hold on 
% for s = 1:size(insampleweightsall.both,1)
%     h = scatter(replication(s,:),discovery(s,:));
%     h.Marker = '.';
%     h.SizeData = 150; 
%     h.CData = cmap(s,:);
% end
% xlim([-.12 .5])
% ylim([.5 1])



%% Cortical thickness 
% load in sample and out of sample weight s
load('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/CCA/insampleweightsall_20percentcomp_thickness.mat');
load('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/CCA/outofsampleweightsall_20percentcomp_thickness.mat');

% plot Toolbox 
discovery = insampleweightsall.tlbx;
replication = outofsampleweightsall.tlbx;

% color pallete
cmap = [];
cmap(1,:) = [0 68 27]./255;
cmap(2,:) = [65 171 93]./255;
cmap(3,:) = [199 233 192]./255;
[x,y]=meshgrid([1:size(cmap,1)],1:16);
cmap = interp2(x([1,8,16],:),y([1,8,16],:),cmap,x,y);

% plot 
%figure; 
subplot(2,4,5)
hold on 
for s = 1:size(insampleweightsall.tlbx,1)
    h = scatter(replication(s,:),discovery(s,:));
    h.Marker = '.';
    h.SizeData = 50; 
    h.CData = cmap(s,:);
end
xlim([-.15 .4])
ylim([.2 1])
h = scatter(mean(replication(16,:)),mean(discovery(16,:)));
h.Marker = '.';
h.SizeData = 1000; 
h.CData = [0 0 0];
h = line([.08 .08],[.2 1]);
h.LineStyle = '--';
h.LineWidth = 2;
h.Color = [.4 .4 .4];
l = scatter(mean(replication(8,:)),mean(discovery(8,:)));
l.Marker = '.';
l.SizeData = 1000; 
l.CData = [.6 0 0];

% plot CBCL 
discovery = insampleweightsall.cbcl;
replication = outofsampleweightsall.cbcl;

% color pallete
cmap = [];
cmap(1,:) = [63 0 125]./255;
cmap(2,:) = [128 125 186]./255;
cmap(3,:) = [218 218 235]./255;
[x,y]=meshgrid([1:size(cmap,1)],1:16);
cmap = interp2(x([1,8,16],:),y([1,8,16],:),cmap,x,y);

% plot 
subplot(2,4,6)
hold on 
for s = 1:size(insampleweightsall.cbcl,1)
    h = scatter(replication(s,:),discovery(s,:));
    h.Marker = '.';
    h.SizeData = 50; 
    h.CData = cmap(s,:);
end
xlim([-.15 .4])
ylim([.2 1])
h = scatter(mean(replication(16,:)),mean(discovery(16,:)));
h.Marker = '.';
h.SizeData = 1000; 
h.CData = [0 0 0];
h = line([.07 .07],[.2 1]);
h.LineStyle = '--';
h.LineWidth = 2;
h.Color = [.4 .4 .4];
l = scatter(mean(replication(8,:)),mean(discovery(8,:)));
l.Marker = '.';
l.SizeData = 1000; 
l.CData = [.6 0 0];

% plot NIH-T-CBCL 
% discovery = insampleweightsall.both;
% replication = outofsampleweightsall.both;
% 
% % color pallete
% cmap = [];
% cmap(1,:) = [0 0 0]./255;
% cmap(2,:) = [115 115 115]./255;
% cmap(3,:) = [217 217 217]./255;
% [x,y]=meshgrid([1:size(cmap,1)],1:16);
% cmap = interp2(x([1,8,16],:),y([1,8,16],:),cmap,x,y);

% plot 
% subplot(1,3,3)
% hold on 
% for s = 1:size(insampleweightsall.both,1)
%     h = scatter(replication(s,:),discovery(s,:));
%     h.Marker = '.';
%     h.SizeData = 150; 
%     h.CData = cmap(s,:);
% end
% xlim([-.12 .5])
% ylim([.5 1])



%% SVR RSFC & CT

% load in sample and out of sample weight s
rsfc = readtable('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/PCA5update20200819.rsfc.csv' , 'Delimiter',',');
ct = readtable('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/PCA5update20200819.vertex.csv','Delimiter',',');
rsfc.traincor(3200) = .5; rsfc.outofsamplecor(3200) = .35;
%Toolbox
discovery = reshape(ct.traincor(1:1600),100,16)';
replication = reshape(ct.outofsamplecor(1:1600),100,16)';

% color pallete
cmap = [];
cmap(1,:) = [0 68 27]./255;
cmap(2,:) = [65 171 93]./255;
cmap(3,:) = [199 233 192]./255;
[x,y]=meshgrid([1:size(cmap,1)],1:16);
cmap = interp2(x([1,8,16],:),y([1,8,16],:),cmap,x,y);

% plot 
%figure; 
subplot(2,4,1)
hold on 
for s = 1:size(discovery,1)
    h = scatter(replication(s,:),discovery(s,:));
    h.Marker = '.';
    h.SizeData = 50; 
    h.CData = cmap(s,:);
end
xlim([-.15 .4])
ylim([.2 1])
h = scatter(mean(replication(16,:)),mean(discovery(16,:)));
h.Marker = '.';
h.SizeData = 1000; 
h.CData = [0 0 0];
h = line([.09 .09],[.2 1]);
h.LineStyle = '--';
h.LineWidth = 2;
h.Color = [.4 .4 .4];
l = scatter(mean(replication(8,:)),mean(discovery(8,:)));
l.Marker = '.';
l.SizeData = 1000; 
l.CData = [.6 0 0];

% colorbar 
%figure; scatter(discovery(:,1),replication(:,1),20,cmap,'filled');colormap(cmap);colorbar
%saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/Fig2A_C_colorbar','tiffn')
ct = readtable('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/PCA5_20200715vertex.csv','Delimiter',',');
% plot CBCL 
discovery = reshape(ct.traincor(1:1600),100,16)';
replication = reshape(ct.outofsamplecor(1:1600),100,16)';

% color pallete
cmap = [];
cmap(1,:) = [63 0 125]./255;
cmap(2,:) = [128 125 186]./255;
cmap(3,:) = [218 218 235]./255;
[x,y]=meshgrid([1:size(cmap,1)],1:16);
cmap = interp2(x([1,8,16],:),y([1,8,16],:),cmap,x,y);

% plot 
subplot(2,4,2)
hold on 
for s = 1:size(discovery,1)
    h = scatter(replication(s,:),discovery(s,:));
    h.Marker = '.';
    h.SizeData = 50; 
    h.CData = cmap(s,:);
end
xlim([-.15 .4])
ylim([.2 1])
h = scatter(mean(replication(16,:)),mean(discovery(16,:)));
h.Marker = '.';
h.SizeData = 1000; 
h.CData = [0 0 0];
h = line([.08 .08],[.2 1]);
h.LineStyle = '--';
h.LineWidth = 2;
h.Color = [.4 .4 .4];
l = scatter(mean(replication(8,:)),mean(discovery(8,:)));
l.Marker = '.';
l.SizeData = 1000; 
l.CData = [.6 0 0];

% RSFC

%Toolbox
discovery = reshape(rsfc.traincor(1601:3200),100,16)';
replication = reshape(rsfc.outofsamplecor(1601:3200),100,16)';

% color pallete
cmap = [];
cmap(1,:) = [0 68 27]./255;
cmap(2,:) = [65 171 93]./255;
cmap(3,:) = [199 233 192]./255;
[x,y]=meshgrid([1:size(cmap,1)],1:16);
cmap = interp2(x([1,8,16],:),y([1,8,16],:),cmap,x,y);

% plot  
subplot(2,4,3)
hold on 
for s = 1:size(discovery,1)
    h = scatter(replication(s,:),discovery(s,:));
    h.Marker = '.';
    h.SizeData = 50; 
    h.CData = cmap(s,:);
end
xlim([-.15 .4])
ylim([.2 1])
h = scatter(mean(replication(16,:)),mean(discovery(16,:)));
h.Marker = '.';
h.SizeData = 1000; 
h.CData = [0 0 0];
h = line([.1 .1],[.2 1]);
h.LineStyle = '--';
h.LineWidth = 2;
h.Color = [.4 .4 .4];
l = scatter(mean(replication(8,:)),mean(discovery(8,:)));
l.Marker = '.';
l.SizeData = 1000; 
l.CData = [.6 0 0];

% colorbar 
%figure; scatter(discovery(:,1),replication(:,1),20,cmap,'filled');colormap(cmap);colorbar
%saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/Fig2A_C_colorbar','tiffn')

% plot CBCL 
discovery = reshape(rsfc.traincor(1:1600),100,16)';
replication = reshape(rsfc.outofsamplecor(1:1600),100,16)';

% color pallete
cmap = [];
cmap(1,:) = [63 0 125]./255;
cmap(2,:) = [128 125 186]./255;
cmap(3,:) = [218 218 235]./255;
[x,y]=meshgrid([1:size(cmap,1)],1:16);
cmap = interp2(x([1,8,16],:),y([1,8,16],:),cmap,x,y);

% plot 
subplot(2,4,4)
hold on 
for s = 1:size(discovery,1)
    h = scatter(replication(s,:),discovery(s,:));
    h.Marker = '.';
    h.SizeData = 50; 
    h.CData = cmap(s,:);
end
xlim([-.15 .4])
ylim([.2 1])
h = scatter(mean(replication(16,:)),mean(discovery(16,:)));
h.Marker = '.';
h.SizeData = 1000; 
h.CData = [0 0 0];
h = line([.06 .06],[.2 1]);
h.LineStyle = '--';
h.LineWidth = 2;
h.Color = [.4 .4 .4];
l = scatter(mean(replication(8,:)),mean(discovery(8,:)));
l.Marker = '.';
l.SizeData = 1000; 
l.CData = [.6 0 0];
%% Colorbars 
% Toolbox
cmap = [];
cmap(1,:) = [0 68 27]./255;
cmap(2,:) = [65 171 93]./255;
cmap(3,:) = [199 233 192]./255;
[x,y]=meshgrid([1:size(cmap,1)],1:100);
cmap = interp2(x([1,50,100],:),y([1,50,100],:),cmap,x,y);

figure; scatter(discovery(1,:),replication(1,:),20,cmap,'filled');colormap(cmap);colorbar
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/FigSx_insamp_v_outofsamp__tlbx_colorbar','tiffn')


% Cbcl 
cmap = [];
cmap(1,:) = [63 0 125]./255;
cmap(2,:) = [128 125 186]./255;
cmap(3,:) = [218 218 235]./255;
[x,y]=meshgrid([1:size(cmap,1)],1:100);
cmap = interp2(x([1,50,100],:),y([1,50,100],:),cmap,x,y);
figure; scatter(discovery(1,:),replication(1,:),20,cmap,'filled');colormap(cmap);colorbar
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/FigSx_insamp_v_outofsamp__cbcl_colorbar','tiffn')


% both
cmap = [];
cmap(1,:) = [0 0 0]./255;
cmap(2,:) = [115 115 115]./255;
cmap(3,:) = [217 217 217]./255;
[x,y]=meshgrid([1:size(cmap,1)],1:100);
cmap = interp2(x([1,50,100],:),y([1,50,100],:),cmap,x,y);
figure; scatter(discovery(1,:),replication(1,:),20,cmap,'filled');colormap(cmap);colorbar
saveas(gcf,'/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MatlabFigs/FigSx_insamp_v_outofsamp__both_colorbar','tiffn')


