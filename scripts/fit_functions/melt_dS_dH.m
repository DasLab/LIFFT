function [f, p1_name, p2_name, variable_parameter_name ] = melt( temperatures, delS, delH );
% Fit to temperature-dependent melts, with entropy change (delS in cal/mol/K) and 
%  enthalpy change (delH in kcal/mol)
%  as parameters.
%
% Equilibrium between two states, K = 
%            exp[ (/R) * ( delS - delH/T) ]
%
% Frac. folded = K/ (1 + K )
%
% Note: assumes no heat capacity difference
%
% (C) R. Das, Stanford University 2008-2016.

p1_name = 'delta-S';
p2_name = 'delta-H';
variable_parameter_name = 'temperature';
R = 0.001986; % kcal/mol/K

% convert celsius to K
K = exp( (1/R) * (delS/1000 - delH./(temperatures + 273.15)));

pred = 1.0 ./ (1.0 + K);
f = [1.0-pred;pred];
