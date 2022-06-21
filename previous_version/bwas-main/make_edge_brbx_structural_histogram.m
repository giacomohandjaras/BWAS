function Edge_brbx_correlations = make_edge_brbx_structural_histogram(allmats_ordered,factors)
% Function to determine edge-level brain behavior correlations. Each
% ROI is extracted for each subject, averaged, and then
% correlated with behavior across subjects. 
%
% Output: Edge_brbx_correlations = roi vector of brain
% behavior correlations 
%
% Input: allmats_ordered = order subject-level correlation matrices (roi x
%              subject) in which ROIs are already sorted by network 
%              (divisions given by 'partitions').
%        factors = subject x behavioral variable matrix 


% number of factors 
varnum = size(factors,2);

% Initialize 
Edge_brbx_correlations = zeros(size(allmats_ordered,2),varnum);

% Loop thru all factors and networks 
for f = 1:size(factors,2)
    for r = 1:size(allmats_ordered,2)
        Edge_brbx_correlations(r,f) = corr(allmats_ordered(:,r),factors(:,f));   
    end  
end

% Vectorize if multiple variables
if size(Edge_brbx_correlations,2) > 1 
    Edge_brbx_correlations = Edge_brbx_correlations(:);
end

end

