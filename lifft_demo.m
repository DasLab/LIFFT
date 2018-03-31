function lifft_demo( demo_name )
%  rundemo( demo_name )
%
%  Demo of LIFFT software for fitting RNA folding transitions 
%   from chemical mapping data
%
%  Possible demo_name's are: 
%
%  'all'           = run all of them
%  'single_ligand' = one ligand binding (FMN binding to a sensor designed on EteRNA) 
%                       (SHAPE mapping data)
%  'double_ligand',
%  'one_two'       = one/two ligand binding (gly. riboswitch binding to 0,1,2 glycines)
%                       (DMS mapping data)
%  'mg'            = Hill fit to a Mg2+ titration (folding of an SRP internal loop)
%  'melt'          = Temperature dependent unfolding of a small AAAA-tetraloop hairpin 
%                       (DMS mapping data)
%  'melt_dS_dH'     = (Not recommended) Same as melt, with dS/dH as
%                        fit parameters insrtead of Tm/dH.
%  'melt_with_linear_baseline'
%                  = Temperature dependent unfolding of a small AAAA-tetraloop hairpin 
%                       (UV absorbance data measured on several cuvettes at different [RNA])

if ( nargin < 1 )
    help( mfilename );
    return;
end
if ~exist( 'saveFNKT_finetitration_051911_conc_area_peak.mat', 'file' ) 
    fprintf( 'Cannot find demos. Make sure to add LIFFT directory and subfolders in "Set Path..." menu!\n')
    return;
end

switch demo_name
    case 'all'
        mg_demo(); pause;
        single_ligand_demo(); pause
        one_two_demo(); pause;
        melt_demo(); pause;
        melt_dS_dH_demo(); pause;
        melt_with_linear_baseline_demo();
    case 'single_ligand'
        single_ligand_demo()
    case {'double_ligand','one_two'}
        one_two_demo()
    case 'mg'
        mg_demo()
    case 'melt'
        melt_demo()
    case 'melt_dS_dH'
        melt_dS_dH_demo()
    case 'melt_with_linear_baseline'
        melt_with_linear_baseline_demo()
    otherwise
        fprintf( 'Unrecognized demo name. For list of demos, type: \nhelp lifft\n' );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mg_demo()
% Adapted from examples/SRPredesign_Mg_binding/SRPredesign_Mgtitration_script_LIFFTtest.m

load Das_969_SRP_DMS_area_norm_conc.mat area_norm conc
which_res = [3:31 33:61]; % there was an additional feature added in the published analysis. This range is what was actually output to RDAT
which_conc = [1:13];  % 14,15,16 show another transition

% Hill coefficients to scan over. In this case, just assumed 1.
n = [1:0.2:4];

% K_d (dissociation constants) to check. Again, units of uM.
K = 10.^[-2:0.05:-0.5];

fit_res  = [17:22 34:39]; % fit to just loop that folds
%fit_res  = [15:23 33:4]; % fit to all residues
plot_res = [19 20 35 36];

do_centralize = 0;
%which_lanes = [1:16]; % wild type
which_lanes = [49:64]; % double mutant (C6G/C22U)
lifft( area_norm(which_res,which_lanes(which_conc)), conc(which_conc), [], K, n, fit_res,'hill',[],plot_res, do_centralize );
fprintf( '\n\nDemo: Hill fit vs. Mg(2+).\nData are for a Rosetta-redesigned SRP domain from Das et al., 2010. https://doi.org/10.1038/nmeth.1433 \n\n\n' )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function single_ligand_demo()
% adapted from
% examples/eterna_fit_FMN_binding/RD031312_EteRNA_MinFMN_titrations8B_script_LIFFTtest.m

load saveRD031312_EteRNA_MinFMN_titrations8B

% go through a few RNA sequences for which we measured binding affinity to flavin mononucleotide (FMN).
i = 1;

% Concentrations over which we titrated FMN. These happen to be in units of micromolar (uM).
conc = [0 0.05 0.075 0.1 0.125 0.15 0.175 0.20 0.25 0.3 0.4 0.5 0.6 0.7 0.8 1.0 1.25 1.5 2.0 3.0 5.0 7.5 10.0 15.0 30.0 50.0 75 150 200 500 1000 5000];

% Hill coefficients to scan over. In this case, just assumed 1.
n = [1];

% K_d (dissociation constants) to check. Again, units of uM.
K = 10.^[-2:0.1:2];

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

fprintf( '\n\nDemo: Langmuir isotherm vs. [ligand].\nData are for an FMN binding design from the EteRNA project. https://doi.org/10.1073/pnas.1313039111\n\n\n' )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function one_two_demo()
% adapted from examples/FN_KT_glycine_riboswitch_two_ligands/FNKT_finetitration_two_ligands_LIFFTtest.m
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

fprintf( '\n\nDemo: F. nucleatum glycine riboswitch vs. glycine.\nData are DMS mapping, fit to "one_two" function (not Hill!). https://doi.org/10.1021/ja2093508\n\n\n' )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function melt_demo()
% adapted from examples/DMS_melt_AAAA/analyze_121814_Matt2_C5_P1_Elim_393027_DMS_script.m
%  Note: fewer input parameters searched, for speed.

% data is saved in here, along with example analysis.
load save_analyze_121814_Matt2_C5_P1_Elim_393027_DMS.mat

Tm = [30:1:80];
delH = [-100:2:-10];

whichres = [-14:-10, 1:12, 23:27] ; % target hairpin, and flanking GAGUA hairpins.
subset = [1:8 10:24];  % leaving out 9th profile due to experimental issue
plotres = [2 4 10]; % some residues in target stem, just for plotting

[ p1_best, p2_best, log_L, C_state, ~, conc_fine, ~, input_data_renorm, pred_fit_fine_renorm ] = lifft( input_data(:,subset), temperatures_good(subset), seqpos_out, Tm, delH, whichres, 'melt', [], plotres,0 );

fprintf( '\n\nDemo: RNA hairpin unfolding vs. temperature.\nData are DMS mapping for an AAAA-tetraloop hairpin (manuscript in prep.).\n\n\n' )



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function melt_dS_dH_demo()
% adapted from examples/DMS_melt_AAAA/analyze_121814_Matt2_C5_P1_Elim_393027_DMS_script.m
%  Note: fewer input parameters searched, for speed.

% data is saved in here, along with example analysis.
load save_analyze_121814_Matt2_C5_P1_Elim_393027_DMS.mat

delS = [-95:0.5:-70];
delH = [-35:0.2:-20];

whichres = [-14:-10, 1:12, 23:27] ; % target hairpin, and flanking GAGUA hairpins.
subset = [1:8 10:24];  % leaving out 9th profile due to experimental issue
plotres = [2 4 10]; % some residues in target stem, just for plotting

[ p1_best, p2_best, log_L, C_state, ~, conc_fine, ~, input_data_renorm, pred_fit_fine_renorm ] = lifft( input_data(:,subset), temperatures_good(subset), seqpos_out, delS, delH, whichres, 'melt_dS_dH', [], plotres,0 );

fprintf( '\n\nDemo (dH vs. dS): RNA hairpin unfolding vs. temperature.\nData are DMS mapping for an AAAA-tetraloop hairpin (manuscript in prep.).\n\n\n' )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function melt_with_linear_baseline_demo()
% adapted from
% examples/optical_melt_AAAA/run_lifft_on_090115_AAAA14_Caco_5_90_optical_melt.m
%  Note: fewer input parameters searched, for speed.

Tm = [30:1:80];
delH = [-100:5:-10];

min_type = 'melt_with_linear_baseline';
min_frac_err = 0.01;

tag = '090115_AAAA14_Caco_5_90/090115_AAAA14_Caco_5_90_exp';

this_dirname = fileparts( which( 'run_lifft_on_090115_AAAA14_Caco_5_90_optical_melt.m' ) );
tag = [this_dirname,'/',tag]
[d,temperatures,filenames] = read_melt_data( tag );
plotres = [2  10]; % some of the melt curves, just for plotting

[ p1_best, p2_best, log_L, C_state, ~, conc_fine, ~, input_data_renorm, pred_fit_fine_renorm  ]  = lifft( d, temperatures, [], Tm, delH,[],min_type,[],plotres,0,[],min_frac_err);

fprintf( '\n\nDemo: RNA hairpin unfolding vs. temperature.\nData are UV absorbance at 260 nm for an AAAA-tetraloop hairpin at different RNA concentrations (manuscript in prep.).\n\n\n' )
