__SBATCH__

. ../../../../axisem3d/__HPC___env.sh

# mesh
cd ../__MESH_DIR__/
cp ../../../../axisem3d/AxiSEM3D_MESHER_ULVZ_build/axisem3d_mesher ./
mpirun -np __N__ axisem3d_mesher

# solve
cd ../__SOLVE_DIR__/
cp ../../../../axisem3d/AxiSEM3D_SOLVER_build/axisem3d ./
mpirun -np __N__ axisem3d
