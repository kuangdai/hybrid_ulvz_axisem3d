#!/bin/bash
#SBATCH --partition=day
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=__N__
#SBATCH --time=__T__
#SBATCH --output=__J__.log
#SBATCH --error=__J__.err
#SBATCH --job-name=__J__
