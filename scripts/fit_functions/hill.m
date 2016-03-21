function [f, p1_name, p2_name ] = hill( conc, K1, n );
%  [f, p1_name, p2_name ] = hill( conc, K1, n );
%
% Hill fit function for LIFFT:
%
% Frac folded = [ conc/K ]^n  / [ 1  +  (conc/K)^n ];
%
% (C) R. Das, Stanford University 2008-2016.

p1_name = 'K';
p2_name = 'nHill';
pred = (conc/K1).^n ./ (1 + (conc/K1).^n );
f = [1-pred; pred ];