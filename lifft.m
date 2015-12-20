function [ p1_best, p2_best, log_L, C_state, input_data_rescale, conc_fine, pred_fit_fine_rescale ] = lifft( input_data, conc, resnum, param1, param2, whichres, fit_type, C_state_in, plot_res, do_centralize, conc_fine )
%  [ p1_best, p2_best, log_L, C_state, input_data_rescale, conc_fine, pred_fit_fine_rescale ] = lifft( input_data, conc, resnum, param1, param2, whichres, fit_type, C_state_in, plot_res, do_centralize, conc_fine );
%
% LIFFT: Likelihood-informed Fits of Footprinting Titraions
%
% Optimizes lane normalizaton and calculates
%  errors at each residue while doing a grid search over midpoints and apparent Hill coefficients.
%
%  Inputs:
%  input_data = matrix of input structure mapping data, approximately normalized (e.g., by mean intensity in each lane)
%                 must have dimensions of (number of residues x number of lanes ).
%  conc       = concentration of chemical in titration (e.g., [adenine], or [Mg2+] )
%  resnum     = your favorite numbering of residues  (if you give empty set [], you'll get 1:number of columns in input_data )
%  K1_conc    = concentrations to search over. (Default: 10.^[-3.0 : 0.1 :3.0]).
%  param2     = apparent Hill coefficients to search over. (Default: param2 = [0.05:0.05:4] )
%  whichres   = [optional!] just look at this subset of input residues. (Give empty set [] to look at all residues. )
%  fit_type   = string specifying functional form to search: "hill", "double_exp", "one_two" ( default is "hill" )
%  C_state_in = a 'target' set of values for footprinting data for each state. If not specified or [], no target set.
%  plot_res   = which residues, if any, to make a 'nice' Hill plot with.  
%  do_centralize = pre-'normalize' the data based on assumption that some residues stay invariant during the titration (default = 1, i.e., true)
%
% (C) Das lab, Stanford University, 2008-2014
% 

% initialization stuff
if size( conc, 2) < size( conc, 1); conc = conc'; end;
if ~isempty( resnum ) & size( resnum, 2) < size( resnum, 1 ); resnum = resnum'; end;
if ( size( input_data, 2) ~= length( conc ) && size( input_data, 1) == length( conc ) ); input_data = input_data'; end;
if isempty( resnum ); resnum = [1:size( input_data, 1) ]; end;
if size( input_data, 2) ~= length( conc )  ;  fprintf( '\nNumber of input_data rows must equal number of values in conc\n' ); return; end;
if size( input_data, 1) ~= length( resnum );  fprintf( '\nNumber of input_data cols must equal number of values in resnum\n' ); return; end;
if ~exist( 'do_centralize'); do_centralize = 1; end;
if ( do_centralize); input_data = centralize( input_data ); end;
if ~exist( 'plot_res' ); plot_res = []; end;

if ~exist( 'param1' ) | isempty( param1 ); param1 = 10.^[-3.0 : 0.1 :3.0]; end
if ~exist( 'param2' ) | isempty( param2 ); param2 = [0.05:0.05:4]; end;
if ~exist( 'resnum' ) | isempty( resnum ); resnum = [1:size( input_data, 1 ) ]; end;
if exist( 'whichres' ) & ~isempty( whichres )
  for k = 1:length(whichres)
    res_to_fit(k) = find( resnum == whichres(k) );
  end
  input_data = input_data( res_to_fit, : );
  resnum = resnum( res_to_fit );
end
if ~exist( 'fit_type' ); fit_type = 'hill'; end;
if ~exist( 'C_state_in' ); C_state_in = []; end;
if ( size( C_state_in, 1 ) == 1 ); C_state_in = input_data( :, C_state_in )'; end;

log_L_best = -99999999999;
p1_best = 999;
p2_best = 999;

% check the evaluation function (default is 'hill' )
[ f, p1_name, p2_name ] = feval( fit_type, conc, param1(1), param2(1) );

tic
% open matlab slaves if they're not open already. Skip this if user does not have parallelization toolbox.
if exist( 'matlabpool');  if matlabpool( 'size' ) == 0 ;   res = findResource; matlabpool( res.ClusterSize ); end; end


% MAIN LOOP! Grid search. Good stuff is inside run_inner_loop
% Will run in parallel ('parfor') if user has parallelization toolbox.
log_L = zeros( length( param1 ), length( param2 ) );
log_L2 = zeros( length( param1 ), length( param2 ) );
sigma_all = zeros( size(input_data,1)+1, length( param1 ), length( param2 ) );
for i = 1: length( param1)
  fprintf(1,'Doing loop %d out of %d...\n', i, length( param1) );
  if exist( 'parfor' )
    parfor j = 1: length( param2) % this could be parallelized for speed..
      [ log_L(i,j), sigma_all(:,i,j) ] = run_inner_loop( param1( i ), param2( j ),...
						  input_data, conc, fit_type, [], C_state_in );
    end
  else
    for j = 1: length( param2) % this could be parallelized for speed..
      [ log_L(i,j), sigma_all(:,i,j) ] = run_inner_loop( param1( i ), param2( j ),...
						  input_data, conc, fit_type, [], C_state_in );
    end
  end
end




t = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Maximum likelihood point.
% maybe this should be replaced by a local minimum finder?
log_L_best = max(max(log_L));
[ind1, ind2] = find(log_L == log_L_best);
p1_best = param1(ind1);
p2_best = param2(ind2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( 'Optimizing further to find local minimum... \n' );
% originally tried fminbound, but that doesn't work for 2D
%p_fminbnd = fminbnd(  'do_new_likelihood_fit_wrapper_for_fminbnd', [param1(ind1-1) param2(ind2-1)], [param1(ind1+1),param2(ind2+1)], [], input_data, conc, fit_type, [], C_state_in );
if length( param2 ) > 1
  [p_fminbnd, minus_log_L_best] = fminsearch(  'do_new_likelihood_fit_wrapper_for_fminbnd', [p1_best, p2_best], [], input_data, conc, fit_type, [], C_state_in );
  p1_best = p_fminbnd(1); p2_best = p_fminbnd(2);
else
  [p_fminbnd, minus_log_L_best] = fminsearch(  'do_new_likelihood_fit_wrapper_for_fminbnd', [p1_best ], [], input_data, conc, fit_type, [], C_state_in, p2_best );
  p1_best = p_fminbnd(1);
end
log_L_best = -minus_log_L_best;

%fprintf(1,'%s %8.4f. %s %8.4f ==>  LogL %8.4f\n', p1_name, p1_best, p2_name, p2_best, log_L_best);

[logLbest, pred_fit, lane_normalization, sigma_at_each_residue, C_state ] =...
    do_new_likelihood_fit( input_data, conc, p1_best, p2_best, fit_type, [],  C_state_in );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Let's try to do a good job of figuring out errors by looking at likelihood near minimum. Where does it drop by 2?
hold on
if ( length( param1 ) > 1 & length( param2 ) > 1 )
 %[p1_low,p1_high] = get_1D_error( param1, max(log_L'),   p1_best, log_L_best );
 %[p2_low,p2_high] = get_1D_error( param2, max(log_L ),   p2_best, log_L_best );
 [p1_low,p1_high] = get_1D_error( param1, max(log_L') );
 [p2_low,p2_high] = get_1D_error( param2, max(log_L ) );
else
 [p1_low,p1_high] = get_1D_error( param1, log_L, p1_best, log_L_best );
 [p2_low,p2_high] = get_1D_error( param2, log_L, p1_best, log_L_best );
end
titlestring = sprintf( '%s = %6.4f + %6.4f - %6.4f\n', p1_name, p1_best, p1_high-p1_best, p1_best - p1_low );
titlestring = [titlestring, sprintf( '%s = %4.2f + %4.2f - %4.2f', p2_name, p2_best, p2_high - p2_best, p2_best - p2_low ) ];
fprintf( [titlestring,'\n'] )

p1_s = [p1_best p1_high-p1_best p1_best - p1_low ];
p2_s = [p2_best p2_high - p2_best p2_best - p2_low];

% Create some smooth fit curves for pretty plots.
log_10_conc = log( conc( find( conc > 0 ) ) ) / log( 10 );
if ~exist( 'conc_fine') | isempty( conc_fine ) 
  if fit_type == 'melt';
    conc_fine = [-10:0.5:110];
  else
    conc_fine = 10.^[(min(log_10_conc)-0.5) : 0.01 : (max(log_10_conc)+0.5) ]; 
  end
end;
f = feval( fit_type, conc_fine, p1_best, p2_best);
pred_fit_fine = C_state'*f;

% Make a pretty plot.
figure(1)
clf
if ( length( param1 ) > 1 & length( param2 ) > 1 )
make_logL_contour_plot( log_L, param1, param2, p1_name, p2_name, p1_best, p2_best );
else
  semilogx( param1, log_L );
end
title( titlestring );

% Make some more pretty plots.
figure(2)
plot_titration_data( input_data, resnum, conc, pred_fit, sigma_at_each_residue, lane_normalization, conc_fine, pred_fit_fine, fit_type );
title( titlestring );

% plot fits and residuals as 'gray plots' too.
figure(4); clf
subplot(1,3,1);
normfactor = mean(mean( input_data ) )/40;
data_lane_norm = input_data*diag(lane_normalization);
image( data_lane_norm/normfactor ); title( 'input data' )
set(gca,'linew',2,'fontsize',14,'fontw','bold','yticklabel',resnum,'ytick',[1:size(input_data,1)]);
set(gca,'linew',2,'fontsize',11,'fontw','bold','xticklabel',conc,'xtick',[1:size(input_data,2)]);
xticklabel_rotate


subplot(1,3,2);
image( pred_fit/normfactor ); title( 'fits' )
set(gca,'linew',2,'fontsize',14,'fontw','bold','yticklabel',resnum,'ytick',[1:size(input_data,1)]);
set(gca,'linew',2,'fontsize',11,'fontw','bold','xticklabel',conc,'xtick',[1:size(input_data,2)]);
xticklabel_rotate

subplot(1,3,3);
image( abs( pred_fit - data_lane_norm)/normfactor ); title( 'abs(residuals)' )
colormap( 1 - gray(100) );
set(gca,'linew',2,'fontsize',14,'fontw','bold','yticklabel',resnum,'ytick',[1:size(input_data,1)]);
set(gca,'linew',2,'fontsize',11,'fontw','bold','xticklabel',conc,'xtick',[1:size(input_data,2)]);
xticklabel_rotate
set(gcf, 'PaperPositionMode','auto','color','white');


% If user has specified 'plot_res' make a plot specifically focused on the data at those residues.
input_data_rescale = []; pred_fit_fine_rescale = [];

if length( plot_res ) > 0; [input_data_rescale, pred_fit_fine_rescale] = make_plot_res_plot( C_state, input_data, lane_normalization, plot_res, conc, resnum, conc_fine, pred_fit_fine,  titlestring, fit_type ); end;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p_low,p_high] = get_1D_error( param, log_L, p_best, log_L_best );

if length( param ) < 2
  p_low  = param;
  p_high = param;
  return;
end

%clf; plot( param, log_L ); pause;

% if user has specified additional point (e.g., best point), stick it into param vector
% this is experimental, and leads to an artefact in 2D scans.
if exist( 'p_best', 'var' ) & exist( 'log_L_best', 'var' )
  insert_pos = 1;
  while param( insert_pos ) < p_best; insert_pos = insert_pos+1; end;
  
  %clf; plot( param, log_L, 'ko-' ); hold on
  param = [ param( [1:(insert_pos-1)] ), p_best    , param( [insert_pos:end]) ];
  log_L = [ log_L( [1:(insert_pos-1)] )', log_L_best, log_L( [insert_pos:end])' ]';

  %plot( param, log_L, 'ro-' ); hold off; pause;
end

[ log_L_max, max_idx ] = max( log_L );
param( max_idx );

log_L_cutoff = log_L_max - 2.0;

p_low = param( max_idx );
if ( max_idx > 1 )
  idx = max_idx-1;
  while (idx > 1 & log_L(idx) < log_L( idx+1) )
    idx = idx - 1;  
  end
  %p_low = param( idx+1 );
  p_low = interp1( log_L( idx: max_idx ), param( [idx : max_idx] ), log_L_cutoff, 'cubic',NaN );
end

p_high = param( max_idx );
if ( max_idx < length( param ) )
  idx = max_idx+1;
  while (idx < length(param) & log_L(idx) < log_L( idx-1 ) )
    idx = idx + 1;  
  end
  %p_high = param( idx-1 );
  p_high = interp1( log_L( max_idx:idx ), param( [ max_idx : idx ] ), log_L_cutoff, 'cubic',NaN );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [logL, sigma_vector] = run_inner_loop( p1, p2, ...
						input_data, conc, fit_type, lane_normalization_in, C_state_in )
    
[logL, pred_fit, lane_normalization, sigma_at_each_residue] = ...
    do_new_likelihood_fit( input_data, conc, p1, p2, fit_type, lane_normalization_in, C_state_in );

sigma_vector = [sigma_at_each_residue' std( lane_normalization )];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [input_data_rescale, pred_fit_fine_rescale] = make_plot_res_plot( C_state, input_data, lane_normalization, plot_res, conc, resnum, conc_fine, pred_fit_fine, titlestring, fit_type ); 

data_renorm = input_data * diag(lane_normalization);

figure(5); clf;
colorcode = jet( length( plot_res ) );
% don't allow pure yellow or light green!
colorcode(:,2) = colorcode(:,2)/2;

for m = 1:length( plot_res )
  i = find( resnum == plot_res( m ) );
  input_data_rescale(:,m) = (data_renorm(i,:) - C_state(1,i))/(C_state(2,i)-C_state(1,i));
  plot( conc, input_data_rescale(:,m), 'o', 'color', colorcode(m,:), 'markerfacecolor',colorcode(m,:) ); hold on;
end

for m = 1:length( plot_res )
  i = find( resnum == plot_res( m ) );
  pred_fit_fine_rescale = (pred_fit_fine(i,:) - C_state(1,i))/(C_state(2,i)-C_state(1,i));
  plot( conc_fine, pred_fit_fine_rescale, '-', 'color', colorcode(m,: ), 'linew',2 );     
end
hold off;
set(gca,'fontweight','bold','fontsize',12,'linew',2);
legend( num2str( plot_res' ), 4 )
if ( fit_type == 'melt' ) 
  variable_parameter_name = 'Temperature'; 
  set(gca,'xscale','lin');
else
  variable_parameter_name = 'Concentration';
  set(gca,'xscale','log');
end;
xlabel( variable_parameter_name ); ylabel( 'Fraction transition' );
xlim( [min( conc_fine ) max( conc_fine ) ] );
set(gcf, 'PaperPositionMode','auto','color','white');
title( titlestring );

