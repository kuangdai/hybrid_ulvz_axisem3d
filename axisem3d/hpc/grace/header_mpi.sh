#!/bin/bash
#SBATCH --partition=day
#SBATCH --nodes=__NODES__
#SBATCH --ntasks-per-node=__N_PER_NODE__
#SBATCH --ntasks=__N__
#SBATCH --time=__T__
#SBATCH --output=__J__.log
#SBATCH --error=__J__.err
#SBATCH --job-name=__J__
