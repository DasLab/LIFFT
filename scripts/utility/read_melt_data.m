function [d, temperatures, filenames] = read_melt_data( tag );
files = dir( [tag,'*.txt'] );
d = []; count = 0;
for i = 1:length( files )
  filename = [ fileparts( tag ), '/', files(i).name ];
  if ~isempty(strfind( filename, '_0.txt' )); continue; end
  if ~isempty(strfind( filename, '_0uM.txt' )); continue; end
  if ~isempty(strfind( filename, 'Buff' )); continue; end
  count = count+1;
  fprintf( 'Reading in %s\n',filename );
  d_file = load( filename );
  temperatures = d_file(:,1);
  d(:,count) = d_file(:,2);
  filenames{count} = filename;
end

