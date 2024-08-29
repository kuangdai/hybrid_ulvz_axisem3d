__SBATCH__

. ../../../../axisem3d/__HPC___env.sh

# solve
cp ../../../../axisem3d/AxiSEM3D_SOLVER_build/axisem3d ./
mpirun -np __N__ axisem3d
