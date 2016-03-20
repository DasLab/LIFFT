% Read in workspace with data.
% These data correspond to two constructs presenting the FMN aptamer in minimal context.
% From:
%  Lee, J., Kladwang, W., Lee, M., Cantu, D., Azizyan, M., Kim, H., Limpaecher, A., Yoon, S., Treuille, A., Das, R., and EteRNA Participants (2014) 
%" RNA design rules from a massive open laboratory" 
%  Proceedings of the National Academy of Sciences U.S.A. 111 (6) : 2122 - 2127 

load saveRD031312_EteRNA_MinFMN_titrations8B

% go through a few RNA sequences for which we measured binding affinity to flavin mononucleotide (FMN).
which_sets = [1 : num_sequences ]; 

% Concentrations over which we titrated FMN. These happen to be in units of micromolar (uM).
conc = [0 0.05 0.075 0.1 0.125 0.15 0.175 0.20 0.25 0.3 0.4 0.5 0.6 0.7 0.8 1.0 1.25 1.5 2.0 3.0 5.0 7.5 10.0 15.0 30.0 50.0 75 150 200 500 1000 5000];

% Hill coefficients to scan over. In this case, just assumed 1.
n = [1];

% K_d (dissociation constants) to check. Again, units of uM.
K = 10.^[-2:0.1:2];

% Simple loop over data sets.
for i = which_sets
  % Number of concentrations (not needed in general, but here I wanted to reverse the order of the input residues, due
  %  to an old convention where the HiTRACE analysis package gave back quantitated chemical mapping data
  %  from 3' to 5'.)
  N = size( area_peak{i},1 );

  % 'Good points' -- over which concentrations were data acquired? (in this case, all 32 conditions were OK). 
  gp = [1:32];

  % Which residues to fit over ('whichres').
  whichres = [15:22];

  plotres = 18;
  lifft( area_peak{i}( :, gp ), conc( gp ), [N: -1 : 1], K, n, whichres, 'hill', [], plotres, 1, [], 0.01 );
   
  pause;
end

% Figure 3: fit curves to signal strength.
% Figure 4: data vs. fits vs. absolute residuals -- should be white (low residuals) on right hand plot.
% Figure 5: The fits currently yield maximum likelihood Kd's of about 0.6 uM and 3.3 uM.

