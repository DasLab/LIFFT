# LIkelihood-based Fits of Folding Transitions (LIFFT)
## What this is
Folding of biomolecules can be induced by cooling temperature, adding small molecule partners, or changing ionic conditions. 

Several techniques now exist to probe these transitions at different places on a molecule. For RNA in particular, chemical mapping
measurements return data at each nucleotide of the molecule, as a function of folding condition: a 2D matrix.

Learning from these data requires modeling them, which in turn requires defining a quantitative model and fitting thermodynamic parameters.
These analyses  have previously required _ad hoc_ assumptions about normalization, which residues to fit, how to estimate errors, etc.

LIFFT defines a likelihood function for these rich data, with thermodynamic relationships vs. temperature, ionic, or ligand condition based on commonly used forms (Hill fits, etc.). 
Errors at each position are inferred from the data themselves, so that positions that vary a lot actually get assigned a high error, while positions that cleanly follow the transition are given lower errors and dominate the fit.
Also, the analysis does a grid search over parameter values and returns contour plots of the likelihood values. So its easy to assess
if parameter uncertainties are correlated to each other, or are highly non-Gaussian.

Written originally by Rhiju Das, Stanford University. original scripts written in 2008, continually expanded to present by Das Lab. 
See publications below for more information on mathematical form.

## Getting started
• Download this package, or clone it.  
• Add the 'scripts' directory to your MATLAB path. Click 'Set Path', and 'Add with Subfolders...', then select 'scripts'.

## Examples
In `examples/eterna_fit_FMN_binding' is an example from the Eterna PNAS 2014 paper, fitting the binding affinity of flavin mononucleotide (FMN) to a couple of  RNA sequences designed to present the FMN aptamer.   

In MATLAB, change working directory to this directory. Then take a look at the example script 'RD031312_EteRNA_MinFMN_titrations8B_script_LIFFTtest.m'. It reads in the
data from a workspace and then runs a loop over two data sets of FMN binding to different RNAs -- run LIFFT on each.
```
open RD031312_EteRNA_MinFMN_titrations8B_script_LIFFTtest
RD031312_EteRNA_MinFMN_titrations8B_script_LIFFTtest
```
Example output, on second data set:

<img width="513" alt="screen shot 2016-03-20 at 11 14 18 am" src="https://cloud.githubusercontent.com/assets/1672331/13906186/bba387ec-ee8d-11e5-8dfc-76b192dfeb75.png">
<img width="528" alt="screen shot 2016-03-20 at 11 14 26 am" src="https://cloud.githubusercontent.com/assets/1672331/13906189/bd92d602-ee8d-11e5-8cc3-55435860bd4a.png">
<img width="494" alt="screen shot 2016-03-20 at 11 14 33 am" src="https://cloud.githubusercontent.com/assets/1672331/13906190/c09f2e22-ee8d-11e5-9245-801aa19e62ce.png">


## More detailed documentation of `lifft` function

To get more detailed documentation, type `help lifft` in MATLAB. Here is a current example:

```
function [ p1_best, p2_best, log_L, C_state, input_data_rescale, conc_fine, pred_fit_fine_rescale, input_data_renorm, pred_fit_fine ] = lifft( input_data, conc, resnum, param1, param2, whichres, fit_type, C_state_in, plot_res, do_centralize, conc_fine, min_frac_error, do_lane_normalization )
%  [ p1_best, p2_best, log_L, C_state, input_data_rescale, conc_fine, pred_fit_fine_rescale, input_data_renorm, pred_fit_fine  ] = lifft( input_data, conc, resnum, param1, param2, whichres, fit_type, C_state_in, plot_res, do_centralize, conc_fine, min_frac_error, do_lane_normalization );
%
% LIFFT: Likelihood-informed Fits of Footprinting Titraions
%
% Optimizes lane normalizaton and calculates
%  errors at each residue while doing a grid search over midpoints and apparent Hill coefficients.
%
%  Inputs:
%
%  input_data = matrix of input structure mapping data, approximately normalized (e.g., by mean intensity in each lane)
%                 must have dimensions of (number of residues x number of lanes ).
%  conc       = concentration of chemical in titration (e.g., [adenine], or [Mg2+] ), or temperature (in C, for melts).
%  resnum     = your favorite numbering of residues  (if you give empty set [], you'll get 1:number of columns in input_data )
%  param1     = parameter 1 values (e.g., concentrations) to search over. (Default: 10.^[-3.0 : 0.1 :3.0]).
%  param2     = parameter 2 values (e.g., apparent Hill coefficients) to search over. (Default: param2 = [0.05:0.05:4] )
%  whichres   = [optional!] just look at this subset of input residues. (Give empty set [] to look at all residues. )
%  fit_type   = string specifying functional form to search: "hill", "double_exp", "one_two", "melt" ( default is "hill" )
%
%  Less commonly used inputs...
%  C_state_in = a 'target' set of values for footprinting data for each state. If not specified or [], no target set.
%  plot_res   = which residues, if any, to make a 'nice' Hill plot with.  
%  do_centralize  = pre-'normalize' the data based on assumption that some residues stay invariant during the titration (default = 1, i.e., true)
%  conc_fine      = finely spaced concentrations to use when plotting. If not specified or [], use default.
%  min_frac_error = minimum assumed relative error in points (default is 0.1)
%  do_lane_normalization = try to fit lane normalization for each parameter value. (default is 1)
%
% Outputs:
%  p1_best = best fit for parameter 1 (e.g., K_d).
%  p2_best = best fit for parameter 2 (e.g., hill coefficient).
%  log_L   = matrix of log-likelihood values over input param1, param2.
%  input_data_rescale = input_data over plot_res, rescaled from 0 to 1 (as shows up in plot_res plot)
%  conc_fine             = concentrations assumed in making a fine curve that goes through input_data_rescale
%  pred_fit_fine_rescale = values from 0 to 1 in fine curve fit through input_data_rescale
%  input_data_renorm     = all input_data, lane normalizations applied
%  pred_fit_fine_renorm  = predicted data, lane normalizations applied
%
```
If you want default values for any inputs, you can typically type in the empty set `[]` for the variable (see above).

## Tips and troubleshooting
• If you see problems with NaN showing up in data (or plots not showing up at all), that's often due to use of thermodynamic parameters 
that lead to undefined 'fraction folded' values, e.g., negative numbers for estimated binding affinities. Adjust the input values
for the thermodynamic parameters.  
• LIFFT, by default, runs a prenormalization ('centralize') on the data based on identifying residues that appear to change the least. If you suspect that this is causing an issue, set the input `do_centralize` to 0.  
• If you want to totally turn off lane normalization during the fit, set `do_lane_normalization` to 0.  
• If you see the fits do not appear to be going through the points cleanly that might be because LIFFT, by default, assumes that the typical error at each residue is > 10%. That level of uncertainty is often the case for chemical mapping, but for optical melts with UV absorbance readout, the error can be much less. In that case, set `min_frac_error` to be lower, e.g. 0.01.  
• As of 2015, LIFFT assumes by default that the data start almost totally in one state and go to another. E.g., almost fully folded to almost fully unfolded. This is assumes as a 'baseline prior' favoring 'fraction folded' near 0 and 1  at beginning and end, respectively. IF you don't want to use this assumption, you'll have to go into the script `do_new_likelihood_fit.m` and find the line `BASELINE_PRIOR = 1`. Change it to `BASELINE_PRIOR = 0`.  
• If you see an error related to `parfor` or similar, that's due to changes in the parallelization toolbox in different MATLAB versions. You can go into the LIFFT function, and comment out the `parfor` block in favor of the `for` block.

## Publications using this package
### For fitting Mg(2+)-dependent folding transitions (incl. Hill fits)
Frederiksen, J.K., Li, N.S., Das, R., Herschlag, D., and Piccirilli, J.A. (2012) "Metal-ion rescue revisited: Biochemical detection of site-bound metal ions important for RNA folding" RNA 18 (6) : 1123 - 1141.
[Paper](https://daslab.stanford.edu/site_data/pub_pdf/2012_Frederiksen_RNA.pdf)

### Standard Hill fits to single ligands
Lee, J., Kladwang, W., Lee, M., Cantu, D., Azizyan, M., Kim, H., Limpaecher, A., Yoon, S., Treuille, A., Das, R., and EteRNA Participants (2014) 
"RNA design rules from a massive open laboratory" 
Proceedings of the National Academy of Sciences U.S.A. 111 (6) : 2122 - 2127.
[Paper](https://daslab.stanford.edu/site_data/pub_pdf/2014_Lee_PNAS.pdf)

### For fitting two-ligand binders (e.g., glycine riboswitch)
Kladwang, W., Chou, F.-C., and Das, R. (2012) "Automated RNA structure prediction uncovers a kink-turn linker in double glycine riboswitches" Journal of the American Chemical Society 134 (3) : 1404 - 1407 
[Paper](https://daslab.stanford.edu/site_data/pub_pdf/2012_Kladwang_JACS.pdf)

### For fitting 'melts' 
Kladwang, W., Seetin, M., Becka, A. Das, R. Standardization of chemical mapping enables accurate quantitation of RNA energetics, in prep.


