__SBATCH__

. ../../../../axisem3d/__HPC___env.sh

# mesh
cd ../@@_mesh_3d/
cp ../../../../axisem3d/AxiSEM3D_MESHER_ULVZ_build/axisem3d_mesher ./
mpirun -np __N__ axisem3d_mesher

# info
cd ../@@_info_3d/
cp ../../../../axisem3d/AxiSEM3D_SOLVER_build/axisem3d ./
mpirun -np __N__ axisem3d

# merge
out=./output
cat $out/WJ_BOX/FLUID_QUAD_KEY_NR.rank*.txt > $out/WJ_BOX/FLUID_QUAD_KEY_NR.txt
cat $out/WJ_BOX/FLUID_REC_KEY_SPZ.rank*.txt > $out/WJ_BOX/FLUID_REC_KEY_SPZ.txt
cat $out/WJ_BOX/SOLID_QUAD_KEY_NR.rank*.txt > $out/WJ_BOX/SOLID_QUAD_KEY_NR.txt
cat $out/WJ_BOX/SOLID_REC_KEY_SPZ.rank*.txt > $out/WJ_BOX/SOLID_REC_KEY_SPZ.txt
rm $out/WJ_BOX/FLUID_QUAD_KEY_NR.rank*.txt
rm $out/WJ_BOX/FLUID_REC_KEY_SPZ.rank*.txt
rm $out/WJ_BOX/SOLID_QUAD_KEY_NR.rank*.txt
rm $out/WJ_BOX/SOLID_REC_KEY_SPZ.rank*.txt
