function [Tm, delH ] = fit_optical_melt( infile )

if exist( infile, 'dir' )
  files = dir( [infile,'/*.txt'] );
  fileID = fopen(strcat(infile, '.txt'), 'w');
  for i = 1:length( files )
    [Tm(i), delH(i) ] = fit_optical_melt( strcat(infile, '/', files(i).name) );
    fprintf(fileID, '%s\t%f\t%f\n', files(i).name, Tm(i), delH(i));
  end
  return;
end

d = load( infile );
[filepath, filename, fileext] = fileparts(infile);

p_in = [50 -20 d(1,2) 0 1.05*d(1,2) 0 ];
p = fminsearch( 'optical_melt', p_in, optimset(), d(:,1), d(:,2 ) );

[~, absorbance_pred] =  optical_melt(p, d(:,1), d(:,2) );

plot( d(:,1), d(:,2), 'bo' ); hold on
plot( d(:,1), absorbance_pred, 'r-' ); hold off

Tm = p(1);
delH = p(2);
title( sprintf( '%s  Tm = %5.1f delH = %6.1f', filename, Tm, delH),'interpreter','none','fontweight','bold' );

set(gcf, 'PaperPositionMode','auto','color','white');
pdfname = [filepath, '/', filename, '.pdf'];
fprintf( ['Exporting pdf: ', pdfname,'\n'] );
export_fig( pdfname );
