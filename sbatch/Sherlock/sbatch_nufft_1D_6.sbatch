#!/bin/bash
#SBATCH --job-name=nufft_1D_6
#SBATCH --output=out/out_nufft_1D_6.out
#SBATCH --error=err/err_nufft_1D_6.err
#SBATCH --time=16:00:00
#SBATCH --partition=bigmem
#SBATCH --qos=bigmem
#SBATCH --nodes=1
#SBATCH --mem=1500000
#################

module load matlab
cd ../../test
matlab -nojvm -r "batch_nufft_1D_6;quit"
