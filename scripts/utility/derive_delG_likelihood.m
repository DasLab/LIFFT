function [ logL_vs_delG, delG_values ] = derive_delG_likelihood( logL, Tm, delH, delG_values );
% derive_delG_likelihood( logL, Tm, delH );
% TODO generalize to any fit type.

% Basic idea:
%  Have logL vs. two 'coordinates', Tm and delH
%  Do a coordinate transformation, switching Tm_fine for delG37
%  Then maximize over remaining delH coordinate.

if ~exist( 'delG_values' ) delG_values = [-10.0:0.05:10.0]; end;
Tm_fine = [min(Tm):0.01:max(Tm)];

%[delHgrid,Tmgrid] = meshgrid( delH, Tm );
%delG37 = delHgrid .* ( (Tmgrid - 37.0 ) ./ (Tmgrid + 273.15));

[delGgrid,Tmgrid] = meshgrid( delG_values, Tm_fine );
delHgrid = delGgrid./( (Tmgrid - 37.0 ) ./ (Tmgrid + 273.15));

logL_transform = interp2( delH, Tm, logL, delHgrid, Tmgrid );
logL_vs_delG = max( logL_transform, [], 1 );