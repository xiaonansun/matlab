#!/bin/bash
#$ -cwd
/opt/hpc/pkg/MATLAB/R2015b/bin/matlab -nodesktop -r "processCaImagingMCPnev('$1');quit"
