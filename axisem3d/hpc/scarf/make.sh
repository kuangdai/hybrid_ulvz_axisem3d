mkdir -p AxiSEM3D_$1_build && cd $_

rm -rf ./* && cmake -DCMAKE_BUILD_TYPE=release \
-DEIGEN3_ROOT=$HOME/eigen-3.4.0 \
-DBOOST_ROOT=/apps20/sw/amd/Boost/1.77.0-GCC-11.2.0 \
-DFFTW_ROOT=/apps20/sw/amd/FFTW/3.3.10-GCC-11.2.0 \
-DMETIS_ROOT=$HOME/miniconda \
-DNETCDF_ROOT=/apps20/sw/amd/netCDF/4.8.1-gompi-2021b \
-DSimpleBinStream_ROOT=$HOME/simplebinstream \
-DCMAKE_C_COMPILER=/apps20/sw/amd/OpenMPI/4.1.1-GCC-11.2.0/bin/mpicc \
-DCMAKE_CXX_COMPILER=/apps20/sw/amd/OpenMPI/4.1.1-GCC-11.2.0/bin/mpicxx \
-DCMAKE_Fortran_COMPILER=/apps20/sw/amd/OpenMPI/4.1.1-GCC-11.2.0/bin/mpif90 \
-DCMAKE_C_FLAGS="-march=znver2" \
-DCMAKE_CXX_FLAGS="-march=znver2" \
../AxiSEM3D_$1/SOLVER

make -j4

cd ..
