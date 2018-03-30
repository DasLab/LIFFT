# LIkelihood-based Fits of Folding Transitions (LIFFT)
## What this is
Folding of biomolecules can be induced by cooling temperature, adding small molecule partners, or changing ionic conditions. 

Several techniques now exist to probe these transitions at different places on a molecule. For RNA in particular, chemical mapping
measurements return data at each nucleotide of the molecule, as a function of folding condition: a 2D matrix. Spectroscopy approaches (NMR, fluorescence, absorbance) can also return data at multiple frequencies for a molecule undergoing a folding transition.

Learning from these data requires modeling them, which in turn requires defining a quantitative model and fitting thermodynamic parameters.
These analyses require assumptions about normalization, which residues to fit, how to estimate errors, etc. which are not always explicitly defined.

LIFFT makes these fitting assumptions explicit through a Bayesian analysis, with thermodynamic relationships vs. temperature, ionic, or ligand condition based on commonly used forms (Hill fits, etc.). 
Errors at each position are inferred from the data themselves, so that positions that vary a lot actually get assigned a high error, while positions that cleanly follow the transition are given lower errors and dominate the fit.
Also, the analysis does a grid search over parameter values and returns contour plots of the likelihood values. So its easy to assess
if parameter uncertainties are correlated to each other, or are highly non-Gaussian.

Written originally by Rhiju Das, Stanford University. Original scripts written in 2008, continually expanded to present by Das Lab. 
See publications below for more information on mathematical form.

## Getting started
• Download this package, or clone it.  
• Add the 'scripts' directory to your MATLAB path. Click 'Set Path', and 'Add with Subfolders...', then select 'LIFFT'.  
• Check out the demos. Run
```
lifft_demo()
```
to see what to run.

## Documentation
Detailed documentation including explicit examples and tips and troubleshooting are being maintained in the [RiboKit LIFFT page](https://ribokit.github.io/LIFFT/docs/).

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
And: Frederiksen et al. as above [Paper](https://daslab.stanford.edu/site_data/pub_pdf/2012_Frederiksen_RNA.pdf).

### For fitting UV melting data 
Vaidyanathan, P., AlSadhan, I., Merriman, D.K., Al-Hashimi, H.M., and Herschlag, D. (2017) "Pseudouridine and N6-methyladenosine modifications weaken PUF protein/RNA interactions" RNA 23: 611-618.




