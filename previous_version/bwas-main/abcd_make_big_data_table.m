% Merge data tables into one big data table 
addpath(genpath('/data/cn/data1/scripts/CIFTI_RELATED/Resources/read_write_cifti/'))
addpath('/data/nil-bluearc/GMT/Scott/MSC_Subcortical/Scripts')
addpath('/data/nil-bluearc/GMT/Scott/ABCD/Scripts/')
addpath('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/ABCD_pitt_washu/Scripts/')
addpath('/data/nil-bluearc/GMT/Scott/scripts/')
addpath('/data/nil-bluearc/GMT/Scott/NLAtoolbox/')

% Read subject list - all the rest subjects we have 
%RestTable = readtable('/data/nil-bluearc/GMT/Scott/ABCD/SubjectLists/total_subject_list_0327.txt');
RestTable = readtable('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/SubjectLists/total_subjectlist_0708.txt');
%added = readtable('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/SubjectLists/subjects_to_untar_0619.txt');

RestTable.Properties.VariableNames{1} = 'Subject';
subjects = RestTable.Subject;

 %% Make parcel correlation matrix
% % define input dir 
%inputdir = '/data/Daenerys/ABCD/data/abcd_collection3165/';
%numpools = 12;
%make_abcd_parcel_corrmats(subjects,numpools,inputdir,'Gordon');
% 
 %% Make vertexwise correlation matrix
% inputdir = '/data/Daenerys/ABCD/data/abcdbids_output/';
% outname = '/data/nil-bluearc/GMT/Scott/ABCD/Vertexcorrmats/Corrmats_Replication';
% numpools = 7;
% make_abcd_vertex_corrmats(subjects,inputdir,outname);

%% Append correlation matrices & motion information  
%load('/data/nil-bluearc/GMT/Scott/ABCD/Parcelcorrmats/SubjectParcelCorrMats.mat')
load('/data/nil-bluearc/GMT/Scott/ABCD/Parcelcorrmats/SubjectParcelCorrMats_withsubcort.mat');
load('/data/nil-bluearc/GMT/Scott/ABCD/FD/SubjectFD.mat')
load('/data/nil-bluearc/GMT/Scott/ABCD/FD/SubjectFD_Incframes.mat')
load('/data/nil-bluearc/GMT/Scott/ABCD/FD/SubjectFrameTotal.mat')
load('/data/nil-bluearc/GMT/Scott/ABCD/FD/Badtmaskidx.mat')
load('/data/nil-bluearc/GMT/Scott/ABCD/Parcelcorrmats/MissingSubIdx.mat')
load('/data/nil-bluearc/GMT/Scott/ABCD/FD/FDtimecourse.mat')
load('/data/nil-bluearc/GMT/Scott/ABCD/Parcelcorrmats/SplitHalfReliability_rest.mat')
load('/data/nil-bluearc/GMT/Scott/ABCD/CorticalThickness/CorticalThickness.mat');
load('/data/nil-bluearc/GMT/Scott/ABCD/CorticalThickness/MissingSubIdx_CortThickness.mat');  
load('/data/nil-bluearc/GMT/Scott/ABCD/CorticalThickness/ParcelThickness.mat');
IncompleteRest = zeros(size(subjects,1),1);

for i = 1:size(RestTable,1)
    Thismat = Corrmat(:,:,i);
    VertexThickness{i,1} = CorticalThickness(:,i);
    ROIThickness{i,1} = ParcelThickness(:,i);
    % check for nans 
    if sum(any(isnan(Thismat))) > 0
        IncompleteRest(i,1) = 1;
    else
        Corrmatcell{i,1} = Thismat;
    end

end

RestTable(:,2) = Corrmatcell;
RestTable(:,3) = num2cell(SubjectFD);
RestTable(:,4) = num2cell(SubjectFD_Incframes);
RestTable(:,5) = num2cell(SubjectFrameTotal);
RestTable(:,6) = num2cell(MissingSubIdx);
RestTable(:,7) = num2cell(Badtmaskidx);
RestTable(:,8) = num2cell(FDtimecourse);
RestTable(:,9) = num2cell(IncompleteRest);
RestTable(:,10) = num2cell(Reliability);
RestTable(:,11) = VertexThickness;
RestTable(:,12) = ROIThickness;
RestTable(:,13) = num2cell(MissingSubIdx_CortThick);
RestTable.Properties.VariableNames{2} = 'Corrmats';
RestTable.Properties.VariableNames{3} = 'FD';
RestTable.Properties.VariableNames{4} = 'FD_Incframes';
RestTable.Properties.VariableNames{5} = 'FrameTotal';
RestTable.Properties.VariableNames{6} = 'MissingRestSubjects';
RestTable.Properties.VariableNames{7} = 'Badtmaskidx';
RestTable.Properties.VariableNames{8} = 'FDtimecourse';
RestTable.Properties.VariableNames{9} = 'IncompleteRest';
RestTable.Properties.VariableNames{10} = 'Reliability';
RestTable.Properties.VariableNames{11} = 'VertexThickness';
RestTable.Properties.VariableNames{12} = 'ROIThickness';
RestTable.Properties.VariableNames{13} = 'MissingThickness';


% Scanner info
%ScannerTable = readtable('/data/nil-bluearc/GMT/Scott/ABCD/SubjectLists/sub_manufac_pairs.txt');
%ScannerSpecTable = readtable('/data/nil-bluearc/GMT/Scott/ABCD/ABCD_ScannerSpecTable.txt');
%ScannerTable = innerjoin(ScannerTable,ScannerSpecTable,'Keys','Scanner_code');
%RestTable = innerjoin(RestTable,ScannerTable,'Keys','Subject');

% Remove subjects with incomplete rest data 
RestTable(IncompleteRest == 1,:) = [];

% Factorization of toolbox
ToolboxFactors = readtable('/data/nil-bluearc/GMT/Scott/ABCD/factorscores_full_bifactor4factors_20190602.txt');
RestTable = innerjoin(RestTable,ToolboxFactors,'Keys','Subject');

% IQ 
IQ = readtable('/data/nil-bluearc/GMT/Scott/ABCD/SubjectLists/ABCD_WISCV.txt');
RestTable = innerjoin(RestTable,IQ,'Keys','Subject');

% % BMI 
% BMI = readtable('/data/nil-bluearc/GMT/Scott/ABCD/bmi_output.txt');
% BMI.x_subjectkey=arrayfun(@(s)(extractAfter(s(1),'NDAR_')),BMI.x_subjectkey);
% BMI.Properties.VariableNames{1} = 'Subject';
% RestTable = innerjoin(RestTable,BMI,'Keys','Subject');
%  
% RT = readtable('/data/nil-bluearc/GMT/Scott/ABCD/abcde/abcd_nback_rt_table.txt');
% RT.x_subjectkey=arrayfun(@(s)(extractAfter(s(1),'NDAR_')),RT.x_subjectkey);
% RT.Properties.VariableNames{1} = 'Subject';
% RestTable = innerjoin(RestTable,RT,'Keys','Subject');

clear Corrmat ParcelThickness CorticalThickness 
%% Behavior - Discovery set 
%CogTable = readtable('/data/nil-bluearc/GMT/Scott/ABCD_Brain_Cog_Paper/Factorscores_Cognition_4026.txt');

% Discovery CBCL & Total Toolbox Subjects
CBCLTable_g1_discovery = readtable('/data/nil-bluearc/GMT/Scott/ABCD/SubjectLists/G1_Discovery/ABCD_2_0_group1_data_cbcl.txt');

% remove subjects with missing data 
Nanidx = zeros(size(CBCLTable_g1_discovery,1),1);
for s = 1:size(CBCLTable_g1_discovery,1)
    Thissub = table2array(CBCLTable_g1_discovery(s,3:end));
    if any(isnan(Thissub))
        Nanidx(s,1) = 1; % Idx to remove subject 
    end
end
disp(['*****    ' num2str(sum(Nanidx)) ' subjects have missing CBCL discovery data and have been removed    *****'])
CBCLTable_g1_discovery(Nanidx==1,:) = [];

%Tlbx
ToolboxTable_g1_discovery = readtable('/data/nil-bluearc/GMT/Scott/ABCD/SubjectLists/G1_Discovery/ABCD_2_0_group1_data_toolboxtotal.txt');
% remove subjects with missing data 
Nanidx = zeros(size(ToolboxTable_g1_discovery,1),1);
for s = 1:size(ToolboxTable_g1_discovery,1)
    Thissub = table2array(ToolboxTable_g1_discovery(s,3:end));
    if any(isnan(Thissub))
        Nanidx(s,1) = 1; % Idx to remove subject 
    end
end
disp(['*****    ' num2str(sum(Nanidx)) ' subjects have missing TLBX discovery data and have been removed    *****'])
ToolboxTable_g1_discovery(Nanidx==1,:) = [];

% Merge behavioral tables
BehavioralTable_g1 = innerjoin(ToolboxTable_g1_discovery,CBCLTable_g1_discovery,'Keys','Subject');
% Merge behavioral table with rest table
MainTable_g1 = innerjoin(RestTable,BehavioralTable_g1,'Keys','Subject');

%% Behavior - Replication set
% CBCL & Total Toolbox Subjects
% CBCL 
CBCLTable_g2_replication = readtable('/data/nil-bluearc/GMT/Scott/ABCD/SubjectLists/G2_Replication/ABCD_2_0_group2_data_cbcl.txt');
% remove subjects with missing data 
Nanidx = zeros(size(CBCLTable_g2_replication,1),1);
for s = 1:size(CBCLTable_g2_replication,1)
    Thissub = table2array(CBCLTable_g2_replication(s,3:end));
    if any(isnan(Thissub))
        Nanidx(s,1) = 1; % Idx to remove subject 
    end
end
disp(['*****    ' num2str(sum(Nanidx)) ' subjects have missing CBCL replication data and have been removed    *****'])
CBCLTable_g2_replication(Nanidx==1,:) = [];

%Tlbx 
ToolboxTable_g2_replication = readtable('/data/nil-bluearc/GMT/Scott/ABCD/SubjectLists/G2_Replication/ABCD_2_0_group2_data_toolboxtotal.txt');
% remove subjects with missing data 
Nanidx = zeros(size(ToolboxTable_g2_replication,1),1);
for s = 1:size(ToolboxTable_g2_replication,1)
    Thissub = table2array(ToolboxTable_g2_replication(s,3:end));
    if any(isnan(Thissub))
        Nanidx(s,1) = 1; % Idx to remove subject 
    end
end
disp(['*****    ' num2str(sum(Nanidx)) ' subjects have missing TLBX replication data and have been removed    *****'])
ToolboxTable_g2_replication(Nanidx==1,:) = [];

% Merge behavioral tables 
BehavioralTable_g2 = innerjoin(ToolboxTable_g2_replication,CBCLTable_g2_replication,'Keys','Subject');
% Merge behavioral table with rest table
MainTable_g2 = innerjoin(RestTable,BehavioralTable_g2,'Keys','Subject');

% Cleanup
clear Toolbox* CBCL* Badtmaskidx BehavioralTable* Corrmatcell IQ Reliability IncompleteRest MissingSubIdx_CortThick ROIThickness VertexThickness inputdir numpools FDtimecourse Nanidx RestTable SubjectF* Thismat Thissub s i MissingSubIdx

%% Save
save('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MainTable_g1_FD.20.mat','MainTable_g1','-v7.3')
save('/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/MainTable_g2_FD.20.mat','MainTable_g2','-v7.3')