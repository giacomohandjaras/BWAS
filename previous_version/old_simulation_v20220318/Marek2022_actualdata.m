clear
clc

% Import Marek 2022 data of Figure 1F (Right): relationship between RSFC
% and psychopathology
marek = readtable('marek2022_figue1F_right.xlsx');

% Convert table to matrix
marek_data = table2array(marek);

% Define the sample sizes
n_obs = [25, 33, 50, 70, 100, 135, 200, 265,...
    375, 525, 725, 1000, 1430, 2000, 2800, 3604];

% Get the number of random resamplings
nboots = size(marek_data,1);

% Get the number of sample sizes
nsamplesizes = size(marek_data,2);

% Number of ROIs analyzed in the Marek et al., 2022 study
n_rois = 333;

% Significance level
alpha_level = 0.05;

% Preallocate array to store uncorrected pvalues of the correlations at
% different sample sizes and for random resamplings
pvals = nan(size(marek_data));

% For each sample size
for s = 1:nsamplesizes
    
    % Determine the number of observations
    N = n_obs(s);
    
    % For each random resampling
    for b = 1:nboots
        
        % Take the abolute value of r. Please note we take the absolute as
        % we will use a parametric method to assess statistical
        % significance, which postulates symmetry of the null distribution.        
        r = abs(marek_data(b,s));
        
        % Convert r value to t
        t = r*sqrt((N-2)/(1-r^2));
        
        % Estimate the two-tailed pvalue of the relationship
        pvals(b,s) = (1 - tcdf(t,(N-2)))*2;
        
    end
    
end

% Find pvalues passing the uncorrected significance level
sign_pvals = pvals<alpha_level;

% Find pvalues passing the Bonferroni corrected significance level
sign_pvalsbonf = pvals<(alpha_level/n_rois);

% Attribute the r sign to pvalues<0.05 uncorrected
sign_pvals_signed = sign(marek_data).*sign_pvals;

% Attribute the r sign to pvalues<0.05 Bonferroni corrected
sign_pvalsbonf_signed = sign(marek_data).*sign_pvalsbonf;

% Count how many uncorrected pvalues are significant for a correlation
% having opposite sign with respect to the true effect
sign_pvals_opposite = sum(sign_pvals_signed<0);

% Count how many uncorrected pvalues are significant for a correlation
% having same sign of the true relationship
sign_pvals_same = sum(sign_pvals_signed>0);

% Catenate same and opposite sign number of pvalues
sign_pvals_all = cat(1,sign_pvals_same,sign_pvals_opposite);

% Count how many bonferroni corrected pvalues are significant for a 
% correlation having opposite sign with respect to the true effect
sign_pvalsbonf_opposite = sum(sign_pvalsbonf_signed<0);

% Count how many bonferroni corrected pvalues are significant for a 
% correlation having same sign of the true relationship
sign_pvalsbonf_same = sum(sign_pvalsbonf_signed>0);

% Catenate same and opposite sign number of pvalues
sign_pvalsbonf_all = cat(1,sign_pvalsbonf_same,sign_pvalsbonf_opposite);

% Create a table of the number of uncorrected p for positive and negative
% correlations
sign_p_table = array2table(sign_pvals_all,...
    'VariableNames',strsplit(num2str(n_obs)),...
    'RowNames',{sprintf('Positive r, p<%.2f (n)',alpha_level),...
    sprintf('Negative r, p<%.2f (n)',alpha_level)});

% Create a table of the number of uncorrected p for positive and negative
% correlations
sign_pcorr_table = array2table(sign_pvalsbonf_all,...
    'VariableNames',strsplit(num2str(n_obs)),...
    'RowNames',{sprintf('Positive r, pBonf<%.2f (n)',alpha_level),...
    sprintf('Negative r, pBonf<%.2f (n)',alpha_level)});

disp(sign_p_table)

disp(sign_pcorr_table)



