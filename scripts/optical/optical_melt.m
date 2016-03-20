function [deviation, absorbance_pred ] = optical_melt( p, temperatures, absorbance );

Tm   = p(1);
delH = p(2);
 a1  = p(3);
da1  = p(4);
 a2  = p(5);
da2  = p(6);

R = 0.001986; % kcal/mol/K

% convert celsius to K
K = exp( (delH/R) * (1/(Tm+273.15) - 1./(temperatures + 273.15)));
pred_folded = 1.0 ./ (1.0 + K);
pred_unfolded   = 1.0 - pred_folded;

absorbance_pred = ( a1 + da1 * temperatures ) .* pred_folded + (a2 + da2*temperatures ) .* pred_unfolded;
deviation = sum( ( absorbance- absorbance_pred ).^2 );