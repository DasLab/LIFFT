function  [ logL, pred_fit, lane_normalization, sigma_at_each_residue, C_state ] = ...
    do_new_likelihood_fit( data, conc, K1, n, fit_type, lane_normalization, C_state_in, beta_C, min_frac_error, baseline_dev )
% [ logL, pred_fit, lane_normalization, sigma_at_each_residue, C_state ] = ...
%    do_new_likelihood_fit( data, conc, K1, n, fit_type, lane_normalization, C_state_in, beta_C  )
%
%  data       = matrix of input structure mapping data, approximately normalized (e.g., by mean intensity in each lane)
%                 must have dimensions of (number of residues x number of lanes ).
%  conc       = concentration of chemical in titration (e.g., [adenine], or [Mg2+] ), or temperature (in C, for melts).
%  K1         = fit parameter 1 (midpoint, in hill fit)
%  n          = fit parameter 2 (apparent Hill coefficient, in hill fit)
%  fit_type   = string specifying functional form to search: "hill", "double_exp", "one_two" ( default is "hill" )
%  lane_normalization = correction factors for each lane (input empty set [] to estimate during likelihood fit)
%  C_state_in = a 'target' set of values for footprinting data for each state. If not specified or [], no target set.
%  beta_C     = parameter governing how close to stay near input C_state_in [not relevant if C_state_in is not inputted!]
%  min_frac_error = minimum assumed relative error in points (default is 0.1)
%  baseline_dev = amount of deviation to allow for start & end fraction folded to deviate from 1 and 0 (default = 0.05)
%
% (C) Das lab, Stanford University, 2008-2016

if ~exist( 'fit_type' ); fit_type = 'hill'; end;
if ~exist( 'beta_C' ) | isempty(beta_C); beta_C = 0.1; end;  % how close to stay near input C_state_in [not relevant if C_state_in is not inputted!]
if ~exist( 'min_frac_error'); min_frac_error = 0.1; end;
if ~exist( 'baseline_dev' ); baseline_dev = 0.05; end;

f = feval( fit_type, conc, K1, n );

numres = size( data, 1 );
numconc = size( data, 2);

numiter = 10; % this seems to converge fast.

if ~exist( 'C_state_in' )
  C_state_in = [];
end

sigma_at_each_residue = ones( 1, numres );

FIT_LANE_NORMALIZATION = 0;
if isempty( lane_normalization )
  FIT_LANE_NORMALIZATION = 1;
  lane_normalization = ones( 1, numconc );
end

SIGMIN_FRAC = min_frac_error;

for n = 1:numiter

  [pred_fit, C_state ] = do_linear_fit_vs_conc( data, f, lane_normalization, C_state_in, sigma_at_each_residue, beta_C );
  
  sigma_at_each_residue = get_sigma_at_each_residue( data, pred_fit, lane_normalization, SIGMIN_FRAC );
  
  if ( FIT_LANE_NORMALIZATION & n < numiter )
    lane_normalization = normalize_lanes( data, pred_fit, sigma_at_each_residue );
  end

end

logL = - numconc * sum( log( sigma_at_each_residue ) );

if exist( 'C_in' ) logL = logL - beta_C * sum( sum( (( C_in - C_state ) * diag( 1./sigma_normalization)).^2 ) ); end;

BASELINE_PRIOR = ( baseline_dev > 0.0 );
if BASELINE_PRIOR
  logL = logL +  0.5 * ( f(1,1)   - 1).^2 / baseline_dev^2;
  logL = logL +  0.5 * ( f(1,end) - 0).^2 / baseline_dev^2;
end

lane_normalization = 1./lane_normalization;

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pred_fit, C ] = do_linear_fit_vs_conc( data, f, lane_normalization, C_in, sigma_at_each_residue, beta );

numres = size( data, 1 );
numconc = size( data, 2);
num_states = size( f, 1 );

data_renorm = data * diag( lane_normalization );

JUST_SCALE_CIN = 0;
ROBUST_FIT = 0;
if exist( 'C_in' ) & ~isempty( C_in ) & JUST_SCALE_CIN

  % this is a special case -- we 'know' what the states look like, and just
   % need to scale them proportionally.
  pred_fit = ( C_in' * f );
  for a = 1:num_states
    pred_fit_contribution = (C_in( a, : )' * f(a, :));
    %clf; plot( pred_fit ); hold on; plot( pred_fit_contribution, 'r'); pause;
    kappa(a) = sum(sum(data_renorm .* pred_fit_contribution)) / sum(sum(pred_fit .* pred_fit_contribution));
  end

  C = diag(kappa) * C_in;

elseif ROBUST_FIT

  C = zeros(num_states, numres);
  for m = 1:numres
    C(:,m) = robustfit(f', data_renorm(m, :), 'fair', 1.345, 'off');
  end

else
  A = f * f';
  % reactivity of each state.
  B =  data_renorm*f';

  if exist( 'C_in' ) & ~isempty( C_in )
    A = A + beta;
    B = B + beta * C_in';
  end
   
  C = A\B';

  % to enforce that C is positive...
  %for m = 1:size( B, 1 );
  %  C(:,m) = lsqnonneg( A, B(m,:)' );
  %end
  
end

pred_fit = C' * f;

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sigma_at_each_residue = get_sigma_at_each_residue( data, pred_fit, lane_normalization, SIGMIN_FRAC );

data_renorm = data * diag( lane_normalization );
deviation_fit = (pred_fit - data_renorm );

sigma_at_each_residue2 = mean( deviation_fit.^2, 2 ) + (SIGMIN_FRAC *mean( data,2)).^2   ;
sigma_at_each_residue = sqrt( sigma_at_each_residue2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ lane_normalization, sigma_at_each_residue] = ...
    normalize_lanes( data, pred_fit, sigma_at_each_residue );

numres = size( data, 1);
numconc = size( data, 2);

% compute the coefficient for the Lagrange multiplier.
for i = 1:numconc
  beta( i ) = 1 / sum( data(:,i) .* pred_fit(:,i) ./ sigma_at_each_residue.^2  );
  data2( i ) = sum( data(:,i) .* data(:,i) ./ sigma_at_each_residue.^2 );
end

lambda = sum( (beta .* data2) - 1 )/sum( beta );

lane_normalization = ( data2 - lambda ) .* beta;


