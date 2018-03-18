function [f, p1_name, p2_name, variable_parameter_name ] = melt( temperatures, Tm, delH );
% Fit to temperature-dependent melts, with melting temperature (Tm) and enthalpy change (delH)
%  as parameters.
%
% Equilibrium between two states, K = 
%            exp[ (delH/R) * ( 1/Tm - 1/T) ]
%
% Frac. folded = K/ (1 + K )
%
% Note: assumes no heat capacity difference
%
% (C) R. Das, Stanford University 2008-2016.

p1_name = 'Tm';
p2_name = 'delta-H';
variable_parameter_name = 'temperature';
R = 0.001986; % kcal/mol/K

% convert celsius to K
K = exp( (delH/R) * (1/(Tm+273.15) - 1./(temperatures + 273.15)));

pred = 1.0 ./ (1.0 + K);
f = [1.0-pred;pred];
