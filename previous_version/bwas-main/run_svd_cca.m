function [outofsampleweights,insampleweights] = run_svd_cca(insamp_brain,insamp_behavior,outofsamp_brain,outofsamp_behavior,binsize,iter)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

outofsampleweights = zeros(length(binsize),iter);
insampleweights = zeros(length(binsize),iter);
subnum = size(insamp_brain,1);

% CCA 
%[~,~,outofsamp_brain] = svd(outofsamp_brain,'econ');
%outofsamp_brain = outofsamp_brain(:,1:numcomp);

% Regress covariates from behavioral data 
% for v = 1:size(outofsamp_behavior,2)
%     [~,~,r] = regress(outofsamp_behavior(:,v),[cov_oos , ones(length(cov_oos),1) ]);
%     outofsamp_behavior(:,v) = r;
% end

% open parallel pool 
z = parpool(length(binsize));

% Loop thru sample sizes 
parfor bin = 1:length(binsize)
    theseweights = zeros(1,iter);
    inweights = zeros(1,iter);
    %if numcomp > 25 
    %    numcomponents = binsize(bin)-10; 
    %else
    %    numcomponents = numcomp;
    %end
    % apply canonical vectors in smaller discovery sample to left out replication sample 
    for j = 1:iter  
        idx = datasample(1:subnum,binsize(bin));
        [coeffs,scores,~,~,explained,mu] = pca(insamp_brain(idx,:));
        ooscorrds = (outofsamp_brain - mu)*coeffs; % Put left out set in same coordinate space 
        % Get 20% of variance explained 
        numcomponents = find(cumsum(explained)>=20);
        numcomponents = numcomponents(1); % Grab index where cumulative var explained reaches 50% 
%         thesebehaviors = [];
%         for v = 1:size(insamp_behavior,2)
%             [~,~,r] = regress(insamp_behavior(idx,v),[cov_insamp(idx) , ones(length(cov_insamp(idx)),1) ]);
%             thesebehaviors(:,v) = r;
%         end
        [A,B,r] = canoncorr(scores(:,1:numcomponents),insamp_behavior(idx,:));
        inweights(1,j) = r(1);
        theseweights(1,j) = corr((ooscorrds(:,1:numcomponents)*A(:,1)),(outofsamp_behavior*B(:,1)));
    end
    insampleweights(bin,:) = inweights;
    outofsampleweights(bin,:) = theseweights;
end
delete(z)
end