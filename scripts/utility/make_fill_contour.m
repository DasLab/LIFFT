function make_fill_contour( logL, Tm, delH, fill_color );
%  make_fill_contour( logL, Tm, delH, fill_color );

[C,h] = contourf( Tm, delH,  logL', [-4 -2]);
set(h,'visible','off')

% first contour
offset = 1;
numvertices = C(2,offset);
gp = find( ~isnan( C(1,offset+[1:numvertices]) ) ) + offset;
h = fill( C(1,gp), C(2,gp), fill_color);
set(h,'edgecolor',fill_color);
set(h,'FaceAlpha',0.0 );

% second contour
offset = offset + numvertices + 1;
if ( offset > size(C,2) ); return; end;
numvertices = C(2,offset);
gp = find( ~isnan( C(1,offset+[1:numvertices]) ) ) + offset;
h = fill( C(1,gp), C(2,gp), fill_color);
set(h,'edgecolor',fill_color);
set(h,'FaceAlpha',0.5 );

