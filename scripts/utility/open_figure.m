function open_figure( name )
h = findobj( 'Type', 'Figure', 'Name', name );
if length(h) == 0
    h = figure('Name',name);
end
figure(h);
set(gcf, 'PaperPositionMode','auto','color','white');
