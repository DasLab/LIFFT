% Read in workspace with data.
% These data correspond to two constructs presenting the FMN aptamer in minimal context.
% From:
%  
% Das, R., Karanicolas, J., and Baker., D. (2010) "Atomic accuracy in predicting and designing noncanonical RNA structure" Nature Methods 7 : 291 - 294
%  doi:10.1038/nmeth.1433
%
% And the data are available at the RNA Mapping Databank at
%  https://rmdb.stanford.edu/detail/SRPDIV_DMS_0001
%

load Das_969_SRP_DMS_area_norm_conc.mat area_norm conc
which_res = [3:31 33:61]; % there was an additional feature added in the published analysis. This range is what was actually output to RDAT
which_conc = [1:13];  % 14,15,16 show another transition

% Hill coefficients to scan over. In this case, just assumed 1.
n = [0.2:0.2:4];

% K_d (dissociation constants) to check. Again, units of uM.
K = 10.^[-2:0.05:-0.5];

fit_res  = [17:22 34:39]; % fit to just loop that folds
%fit_res  = [15:23 33:4]; % fit to all residues
plot_res = [19 20 35 36];


do_centralize = 0;
which_lanes = [49:64]; % double mutant (C6G/C22U)
lifft( area_norm(which_res,which_lanes(which_conc)), conc(which_conc), [], K, n, fit_res,'hill',[],plot_res );
pause

which_lanes = [1:16]; % wild type
lifft( area_norm(which_res,which_lanes(which_conc)), conc(which_conc), [], K, n, fit_res,'hill',[],plot_res );


%%%%%%%%%%%%%%