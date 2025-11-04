mkdir -p AxiSEM3D_$1_build && cd $_

module purge
#module load CMake/3.24.3-GCCcore-12.2.0 Eigen/3.4.0-GCCcore-12.2.0 Boost/1.81.0-iompi-2022b FFTW/3.3.10-iomkl-2022b METIS/5.1.0-GCCcore-12.2.0-32bit netCDF/4.9.0-iompi-2022b
module load CMake/3.24.3-GCCcore-12.2.0 Eigen/3.4.0-GCCcore-12.2.0 Boost/1.81.0-GCC-12.2.0 FFTW/3.3.10-GCC-12.2.0 METIS/5.1.0-GCCcore-12.2.0-32bit netCDF/4.9.0-gompi-2022b
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

