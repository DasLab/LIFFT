function [minus_logL, pred_fit, lane_normalization, sigma_at_each_residue] = ...
    do_new_likelihood_fit_wrapper_for_fminbnd( p, input_data, conc, fit_type, lane_normalization_in, C_state_in, p2_best, min_frac_error, baseline_dev );
% [minus_logL, pred_fit, lane_normalization, sigma_at_each_residue] = ...
%    do_new_likelihood_fit_wrapper_for_fminbnd( p, input_data, conc, fit_type, lane_normalization_in, C_state_in, p2_best, min_frac_error );
%
% Simple wrapper used for optimization within lifft.m
%
% See 'help do_new_likelihood_fit' for more info on inputs.
%
p1 = p(1); 
if length( p ) > 1
  p2 = p(2);
else
  p2 = p2_best;
end
[logL, pred_fit, lane_normalization, sigma_at_each_residue ] = do_new_likelihood_fit( input_data, conc, p1, p2, fit_type, lane_normalization_in, C_state_in, [], min_frac_error, baseline_dev );
minus_logL = -logL;
