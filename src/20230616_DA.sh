#!/bin/bash
#SBATCH --account=def-nameprof-ab
#SBATCH --mem-per-cpu=4G
#SBATCH --time=0-0:10
module load gcc/9.3.0 r/4.2.2
Rscript 20230616_DA.R
