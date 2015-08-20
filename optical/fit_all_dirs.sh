#!/bin/bash

source ~/.bashrc
for f in ${1}/*; do
    /biox3/software/non-free/MATLAB-R2013b/bin/matlab -nodesktop -nosplash -r "fit_optical_melt('$f');exit;"
done
