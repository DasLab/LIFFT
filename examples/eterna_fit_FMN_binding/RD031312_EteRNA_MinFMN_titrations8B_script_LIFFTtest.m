load saveRD031312_EteRNA_MinFMN_titrations8B

which_sets = [1 : num_sequences ]; 

conc = [0 0.05 0.075 0.1 0.125 0.15 0.175 0.20 0.25 0.3 0.4 0.5 0.6 0.7 0.8 1.0 1.25 1.5 2.0 3.0 5.0 7.5 10.0 15.0 30.0 50.0 75 150 200 500 1000 5000];
n = [1];
K = 10.^[-2:0.1:2];

for i = which_sets
  N = size( area_peak{i},1 );
  gp = [1:32];
  resnum = [15:22];
  lifft( area_peak{i}( :, gp ), conc( gp ), [N: -1 : 1], K, n, resnum, 'hill', [],  [18] );
  pause;
end

