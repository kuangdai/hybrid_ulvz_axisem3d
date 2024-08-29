__SBATCH__

. ../../../../axisem3d/__HPC___env.sh

# compute stf
cp ../../../../axisem3d/AxiSEM3D_STF_build/axisem3d_stf ./
./axisem3d_stf __E_LAT__ __E_LON__ __U_LAT__ __U_LON__ FLUID 0 1000000000 __NR_EVT__ __NR_ULVZ__ -10000 10000 ../@@_stations ../@@_solve_prem ../@@_info_3d
./axisem3d_stf __E_LAT__ __E_LON__ __U_LAT__ __U_LON__ SOLID 0 1000000000 __NR_EVT__ __NR_ULVZ__ -10000 10000 ../@@_stations ../@@_solve_prem ../@@_info_3d
