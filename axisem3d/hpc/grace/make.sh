mkdir -p AxiSEM3D_$1_build && cd $_

rm -rf ./* && cmake -DCMAKE_BUILD_TYPE=release \
-DEIGEN3_ROOT=/gpfs/loomis/apps/avx/software/Eigen/3.4.0-GCCcore-10.2.0/include \
-DBOOST_ROOT=/gpfs/loomis/apps/avx/software/Boost/1.74.0-GCCcore-10.2.0 \
-DFFTW_ROOT=/gpfs/loomis/apps/avx/software/FFTW/3.3.8-gompi-2018b \
-DMETIS_ROOT=/gpfs/loomis/apps/avx/software/METIS/5.1.0-GCCcore-7.3.0-32bi \
-DNETCDF_ROOT=/gpfs/loomis/apps/avx/software/netCDF/4.6.1-foss-2018b \
-DSimpleBinStream_ROOT=$HOME/simplebinstream \
../AxiSEM3D_$1/SOLVER

make

cd ..
