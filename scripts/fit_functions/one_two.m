function [f, p1_name, p2_name, p1_unit_name, p2_unit_name, variable_name, variable_unit_name, variable_scale ] = one_two( conc, K1, K2 );
%  [f, p1_name, p2_name ] = one_two( conc, K1, K2 );
%
% Binding isotherm for a molecule that binds up to 2 ligands.
%
%  X <--> X*L <--> X*2L
%     K1       K2      
%
% K1 and K2 are dissociation constants for first and second ligand.
%
% f0 = 1         / [1 + conc/K1 + (conc/K1)(conc/K2)]
% f1 = (conc/K1) / [1 + conc/K1 + (conc/K1)(conc/K2)]
% f2 = (conc/K2) / [1 + conc/K1 + (conc/K1)(conc/K2)]
%
% (C) R. Das, Stanford University 2011-2016.
%

p1_name = 'K_1';
p2_name = 'K_2';
p1_unit_name = ''; % don't really know if its mM or what.
p2_unit_name = '';
variable_name = 'Concentration';
variable_unit_name = '';
variable_scale = 'log';

f = [ ones(1,length(conc));   (conc/K1);   (conc/K2).*(conc/K1) ];
[dummy, Z]  = meshgrid( ones(1,3), sum( f ) );
f = f./Z';