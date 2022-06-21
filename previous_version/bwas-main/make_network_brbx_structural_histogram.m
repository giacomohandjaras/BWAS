function Network_brbx_correlations = make_network_brbx_structural_histogram(allmats_ordered,factors,partitions)
% Function to determine network-level brain behavior correlations. Each
% block (network pair) is extracted for each subject, averaged, and then
% correlated with behavior across subjects. 
%
% Output: Network_brbx_correlations = network x network matrix of brain
% behavior correlations (diagonal represents within network relationships,
% off diagonal upper triangle represented between network relationships 
%
% Input: allmats_ordered = order subject-level correlation matrices (roi x
%              subject) in which ROIs are already sorted by network 
%              (divisions given by 'partitions').
%        factors = subject x behavioral variable matrix 
%        partitions = division of networks 

% number of factors 
varnum = size(factors,2);

% Create new partitions that include a 1 and last ROI (e.g., 333). 
if partitions(1) == 1
        error('Remove the 1 in your partition index!')
end

temp = ones(length(partitions)+1,1);
temp(2:end) = partitions;
partitions = temp;

Network_brbx_correlations = zeros(length(partitions),varnum);

% Loop thru all factors and networks 
for f = 1:size(factors,2)
    for x = 1:length(partitions)
        if x < length(partitions)
            Network_brbx_correlations(x,f) = corr(mean(allmats_ordered(:,partitions(x):partitions(x+1)-1),2),factors(:,f));
        else
            Network_brbx_correlations(x,f) = corr(mean(allmats_ordered(:,partitions(x):333,:),2),factors(:,f));
        end
    end  
end

% Vectorize if multiple variables
if size(Network_brbx_correlations,2) > 1 
    Network_brbx_correlations = Network_brbx_correlations(:);
end

end

