% fine bins
Tm = [30:0.5:80];
delH = [-100:5:-10];

min_type = 'melt_with_linear_baseline';
min_frac_err = 0.01;
do_centralize = 0;
tag = '090115_AAAA14_Caco_5_90/090115_AAAA14_Caco_5_90_exp';

this_dirname = fileparts( which( 'run_lifft_on_090115_AAAA14_Caco_5_90_optical_melt.m' ) );
tag = [this_dirname,'/',tag]
[d,temperatures,filenames] = read_melt_data( tag );
[ p1_best, p2_best, log_L, C_state, ~, conc_fine, ~, input_data_renorm, pred_fit_fine_renorm  ]  = lifft( d, temperatures, [], Tm, delH,[],min_type,[],[1:size(d,2)],do_centralize,[],min_frac_err);

% at end of run, output originally read:
%
% Tm = 54.7899 + 2.0756 - 1.9008
% delta-H = -38.19 + 7.21 - 10.14
%
% After some small changes to default parameters (do_centralize=0)
% 
% Tm = 55.6577 + 2.7891 - 2.3142
% delta-H = -32.35 + 7.67 - 10.99
%

% saving example output in here.
% save save090115_AAAA14_Caco_5_90_exp.mat;

