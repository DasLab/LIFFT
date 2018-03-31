function [f, p1_name, p2_name, p1_unit_name, p2_unit_name, variable_name, variable_unit_name, variable_scale ]  = melt_with_linear_baseline( temperatures, Tm, delH );
% Fit to temperature-dependent melts, with melting temperature (Tm) and enthalpy change (delH)
%  as parameters.
%
% Equilibrium between two states, K = 
%            exp[ (delH/R) * ( 1/Tm - 1/T) ]
%
% Frac. folded, f = K/ (1 + K )
%
% States have populations f and 1-f.
%
% Also outputs two other 'states' with fraction folded T*f and T*(1-f), which
%  tricks LIFFT into also fitting a linear baseline for the folded and unfolded states.
%
% Note: assumes no heat capacity difference
%
% (C) R. Das, Stanford University 2008-2016.

p1_name = 'T_m';
p2_name = '\DeltaH';
p1_unit_name = '{\circ}C';
p2_unit_name = 'kcal/mol';
variable_name = 'Temperature';
variable_unit_name = '{\circ}C';
variable_scale = 'linear';
R = 0.001986; % kcal/mol/K

% convert celsius to K
K = exp( (delH/R) * (1/(Tm+273.15) - 1./(temperatures + 273.15)));

pred = 1.0 ./ (1.0 + K);
f = [1.0-pred; pred];

f = [f; f.*repmat( temperatures, size(f,1), 1 )];
