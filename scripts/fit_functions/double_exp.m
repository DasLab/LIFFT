function [f, p1_name, p2_name, variable_name, variable_unit_name, variable_scale ] = double_exp( times, tau1, tau2 );
%  [f, p1_name, p2_name ] = double_exp( times, tau1, tau2 );
%
% Double exponential fit function for LIFFT
%
% Three states:
%
% Frac folded = 1, exp( - t/tau1), exp( -t/tau2).
%
% Note: not well tested.
%

p1_name = '\tau_1';
p2_name = '\tau_2';
p1_unit_name = 'sec';
p2_unit_name = 'sec';
variable_name = 'Time';
variable_unit_name = 'sec';
variable_scale = 'linear';

f = [ ones(1,length(times)); exp( -times./tau1); exp( -times./tau2 ) ];