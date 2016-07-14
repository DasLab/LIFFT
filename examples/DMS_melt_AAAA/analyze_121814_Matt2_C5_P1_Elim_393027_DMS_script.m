% data is saved in here, along with example analysis.
load save_analyze_121814_Matt2_C5_P1_Elim_393027_DMS.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tm = [30:0.1:80];
delH = [-100:0.5:-10];

whichres = [-14:-10, 1:12, 23:27] ; % target hairpin, and flanking GAGUA hairpins.
subset = [1:8 10:24];  % leaving out 9th profile due to experimental issue
plotres = [2 4 10]; % some residues in target stem, just for plotting

[ p1_best, p2_best, log_L, C_state, ~, conc_fine, ~, input_data_renorm, pred_fit_fine_renorm ] = lifft( input_data(:,subset), temperatures_good(subset), seqpos_out, Tm, delH, whichres, 'melt', [], plotres,0 );

