function [f, p1_name, p2_name ] = double_exp( times, tau1, tau2 );
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

p1_name = 'tau1';
p2_name = 'tau2';

f = [ ones(1,length(times)); exp( -times./tau1); exp( -times./tau2 ) ];