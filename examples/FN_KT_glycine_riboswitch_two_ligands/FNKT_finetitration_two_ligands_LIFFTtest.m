% These data correspond to glycine titrations of the F. nucleatum glycine riboswitch with a
%  kink-turn linker. 
%
% From:
%  
% Kladwang, W., Chou, F.-C., and Das, R. (2012) "Automated RNA structure prediction uncovers a kink-turn linker in double glycine riboswitches" Journal of the American Chemical Society 134 (3) : 1404 - 1407
%  http://pubs.acs.org/doi/abs/10.1021/ja2093508
%
% And the data are available at the RNA Mapping Databank at
%  https://rmdb.stanford.edu/detail/GLYCFN_KNK_0002
%

load saveFNKT_finetitration_051911_conc_area_peak.mat area_peak conc seqpos

whichres = [ -7:158 ];
subset = [1:2 4:32];
conc_fit = conc(subset); conc_to_plot = conc; conc_to_plot(1) = 1e-2;
whichcols = [1:32 ]; whichcols = whichcols( subset );
input_data = area_peak(:,whichcols);
fit_type = 'one_two';
p1 = 10.^[0:0.05:2]; p2 =10.^[2:0.1:6];
plot_res = [43 52 63 135 137 153 ];
lifft( input_data, conc_fit, seqpos,p1, p2,whichres,fit_type,[], plot_res);