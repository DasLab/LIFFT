function [ p1_best, p2_best, log_L, C_state, input_data_rescale, conc_fine, pred_fit_fine_rescale ] = analyze_titration_with_likelihood( input_data, conc, resnum, param1, param2, whichres, fit_type, C_state_in, plot_res, do_centralize, conc_fine )
% Use LIFFT instead!!!

if ~exist( 'conc_fine', 'var' ) conc_fine = []; end;

fprintf( '\nanalyze_titration_with_likelihood.m is deprecated! Use LIFFT instead!\n\n' );
[ p1_best, p2_best, log_L, C_state, input_data_rescale, conc_fine, pred_fit_fine_rescale ] = lifft( input_data, conc, resnum, param1, param2, whichres, fit_type, C_state_in, plot_res, do_centralize, conc_fine );

fprintf( '\nanalyze_titration_with_likelihood.m is deprecated! Use LIFFT instead!\n\n' );
