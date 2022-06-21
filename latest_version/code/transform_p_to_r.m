function r_crit = transform_p_to_r(two_tailed_p,n_comparisons,n_observations)
% coded by Luca Cecchetti & Giacomo Handjaras, IMT-Lucca, Italy
% vers 20220619

%t = r*sqrt((n-2)/(1-r^2))
%p = 2*tcdf(t,(n-2));

% Apply correction for multiple comparisons
p = two_tailed_p/n_comparisons;

% Number of samples for which the critical r is computed
n_samples = numel(n_observations);

% Preallocate an array to store critical r values
r_crit = nan(n_samples,1);

% For each sample size
for i = 1:n_samples
    
    % Get the number of observations
    sample_size = n_observations(i);
    
    % Back-transform pvalue into critical t statistic
    t = abs(tinv(p/2,(sample_size-2)));
    
    % Transform critical t statistic into critical r
    r_crit(i) = sqrt((t^2) / ((t^2) + (sample_size-2)));    
    
end


end



















