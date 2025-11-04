mkdir -p AxiSEM3D_$1_build && cd $_
module load cmake eigen boost cray-fftw metis netcdf4
module list

export CC=mpicc
export CXX=mpicxx
export FC=mpifort
export METIS_ROOT=/mnt/lustre/a2fs-work4/work/y07/shared/libs/core/metis/5.1.0/CRAYCLANG/15.0
export METIS_INCLUDE=$METIS_ROOT/include
export METIS_LIB=$METIS_ROOT/lib
rm -rf ./* && cmake -DCMAKE_BUILD_TYPE=release \
-DEIGEN3_ROOT=/work/y07/shared/libs/core/eigen/3.4.0 \
-DBOOST_ROOT=/mnt/lustre/a2fs-work4/work/y07/shared/archer2-lmod/libs/core/boost/1.81.0 \
-DFFTW_ROOT=/opt/cray/pe/lmod/modulefiles/cpu/x86-rome/1.0/cray-fftw/3.3.10.5 \
-DEIGEN3_INCLUDE_DIR=/work/y07/shared/libs/core/eigen/3.4.0/include \
-DNETCDF_ROOT=/work/y07/shared/python/core/netcdf4/1.6.4 \
-DSimpleBinStream_ROOT=/work/n03/n03/yhdai/axisem3d_hybird/axisem3d/hpc/Archer2/depends/simplebinstream \
-DMETIS_INCLUDE_DIRS=$METIS_INCLUDE \
-DMETIS_LIBDIR=$METIS_LIB \
-DMETIS_LIBRARY_DIRS=$METIS_LIB \
-DMETIS_LIBRARY=$METIS_LIB/libmetis_crayclang.a \
-DMETIS_LIBRARIES=$METIS_LIB/libmetis_crayclang.a \
-DMETIS_WORKS=TRUE \
../AxiSEM3D_$1/SOLVER

make -j ${SLURM_CPUS_ON_NODE}

cd ..

