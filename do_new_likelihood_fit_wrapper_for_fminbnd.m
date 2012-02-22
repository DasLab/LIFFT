function [minus_logL, pred_fit, lane_normalization, sigma_at_each_residue] = ...
    do_new_likelihood_fit_wrapper_for_fminbnd( p, input_data, conc, fit_type, lane_normalization_in, C_state_in, p2_best );

p1 = p(1); 
if length( p ) > 1
  p2 = p(2);
else
  p2 = p2_best;
end
[logL, pred_fit, lane_normalization, sigma_at_each_residue ] = do_new_likelihood_fit( input_data, conc, p1, p2, fit_type, lane_normalization_in, C_state_in );
minus_logL = -logL;
