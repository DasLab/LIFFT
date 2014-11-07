function [f, p1_name, p2_name ] = melt( temperatures, Tm, delH );
p1_name = 'Tm';
p2_name = 'delta-H';
R = 0.001986; % kcal/mol/K

% convert celsius to K
K = exp( (delH/R) * (1/(Tm+273.15) - 1./(temperatures + 273.15)));

pred = 1.0 ./ (1.0 + K);
f = [pred; 1.0-pred];
