function [f, p1_name, p2_name, p1_unit_name, p2_unit_name, variable_name, variable_unit_name, variable_scale ] = hill( conc, K1, n );
%  [f, p1_name, p2_name ] = hill( conc, K1, n );
%
% Hill fit function for LIFFT:
%
% Frac folded = [ conc/K ]^n  / [ 1  +  (conc/K)^n ];
%
% (C) R. Das, Stanford University 2008-2016.

p1_name = 'K';
p2_name = 'n_{Hill}';
p1_unit_name = ''; % don't really know if its mM or what.
p2_unit_name = '';
variable_name = 'Concentration';
variable_unit_name = '';
variable_scale = 'log';

pred = (conc/K1).^n ./ (1 + (conc/K1).^n );
f = [1-pred; pred ];