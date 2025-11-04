mkdir -p AxiSEM3D_$1_build && cd $_

rm -rf ./* && cmake -DCMAKE_BUILD_TYPE=release \
-DEIGEN3_ROOT=/gpfs/loomis/apps/avx/software/Eigen/Eigen/3.4.0-GCCcore-12.2.0 \
-DBOOST_ROOT=/gpfs/loomis/apps/avx/software/Boost/1.81.0-iompi-2022b \
-DFFTW_ROOT=/gpfs/loomis/apps/avx/software/FFTW/FFTW/3.3.10-iomkl-2022b \
-DMETIS_ROOT=/gpfs/loomis/apps/avx/software/GCCcore-12.2.0 METIS/5.1.0-GCCcore-12.2.0 \
-DNETCDF_ROOT=/gpfs/loomis/apps/avx/software/netCDF/4.9.0-iompi-2022b \
-DSimpleBinStream_ROOT=$HOME/simplebinstream \
../AxiSEM3D_$1/SOLVER

make

cd ..
