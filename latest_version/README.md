%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulation included in the latest version of the manuscript (June 2022)
%%% Descriptions of the files/directories included
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%
external_code -> external libraries (e.g., fdr_bh.m, export_fig, etc...) used to perform analyses

%%%%
code -> our set of functions to simulate data, calculate power, type M and S errors, replication rate

%%%%
derivatives -> the folder with figures and workspaces generated by the parent scripts

%%%%
raw_data -> raw materials, downloaded from source figures of the original M&TC paper (https://www.nature.com/articles/s41586-022-04492-9)

%%%%
manuscript -> current version of our manuscript

%%%%
simulation_parent_script.m -> code to perform statistical measures between rsFC components and cognition, as well as ROI-based cortical thickness values and cognitive ability scores (Figure 1B-C-D)

%%%%
hcp_parent_script.m.m -> code to perform statistical measures between cortical thickness and fluid intelligence with HCP data (Figure 1E)

%%%%
power_decrease_parent_script.m -> code to prove that power estimates for the average method depend on the size of the full sample (Figure 2A)

%%%%
inflation_effect_parent_script.m -> code to assess effect size inflation across multiple statistical thresholds (Figure 2B)

%%%%
sample_adjustment_parent_script.m -> code to demonstrate that effect sizes reported in the literature align with true effects measured in the population only when the sample size is large enough (Figure 2C)