% fine bins
Tm = [30:0.1:80];
delH = [-100:5:-10];

min_type = 'melt_with_linear_baseline';
min_frac_err = 0.01;

tag = '090115_AAAA14_Caco_5_90/090115_AAAA14_Caco_5_90_exp';
[d,temperatures,filenames] = read_melt_data( tag );
[ p1_best, p2_best, log_L, C_state, ~, conc_fine, ~, input_data_renorm, pred_fit_fine_renorm  ]  = lifft( d, temperatures, [], Tm, delH,[],min_type,[],[1:size(d,2)]',0,[],min_frac_err);

% at end of run, output should read:
%
% Tm = 54.7899 + 2.0756 - 1.9008
% delta-H = -38.19 + 7.21 - 10.14
%

% saving example output in here.
% save save090115_AAAA14_Caco_5_90_exp.mat;

