mkdir -p AxiSEM3D_$1_build && cd $_
./env.sh
module list

export FFTW_ROOT=${EBROOTFFTW}
export CC=mpicc
export CXX=mpicxx
export FC=mpifort

rm -rf ./* && cmake -DCMAKE_BUILD_TYPE=release \
-DEIGEN3_ROOT=${EBROOTEIGEN} \
-DBOOST_ROOT=${EBROOTBOOST} \
-DFFTW_ROOT=${EBROOTFFTW} \
-DMETIS_ROOT=${EBROOTMETIS} \
-DNETCDF_ROOT=${EBROOTNETCDF} \
-DSimpleBinStream_ROOT=$HOME/simplebinstream \
../AxiSEM3D_$1/SOLVER

make -j ${SLURM_CPUS_ON_NODE}

cd ..

