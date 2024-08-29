mkdir -p AxiSEM3D_$1_build && cd $_

rm -rf ./* && cmake -DCMAKE_BUILD_TYPE=release \
-DEIGEN3_ROOT=/Users/kuangdai/axisem3d_depends/eigen3-dev \
-DBOOST_ROOT=/Users/kuangdai/axisem3d_depends/boost/boost \
-DFFTW_ROOT=/Users/kuangdai/anaconda3/lib \
-DMETIS_ROOT=/Users/kuangdai/anaconda3/lib \
-DNETCDF_ROOT=/Users/kuangdai/anaconda3/lib \
-DSimpleBinStream_ROOT=/Users/kuangdai/axisem3d_depends/simplebinstream \
-DNO_OMP=1 \
../AxiSEM3D_$1/SOLVER

make -j4

cd ..
