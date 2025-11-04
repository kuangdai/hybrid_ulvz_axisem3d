//
//  mpi.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/6/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  MPI interfaces

#ifndef mpi_hpp
#define mpi_hpp

#include <complex>
#include <vector>
#include <map>

#ifndef _SERIAL_BUILD
// mpi.h has extern "C" {} internally
#include <mpi.h>
#else
#define MPI_Request int
#define MPI_REQUEST_NULL 0
#endif

namespace mpi {
    ////////////////////////////// internal //////////////////////////////
    namespace internal {
#ifndef _SERIAL_BUILD
        // group superior
        extern MPI_Comm iCommSuper;
        // group inferior
        extern MPI_Comm iCommInfer;
        // current MPI_Comm
        extern MPI_Comm iCommCurrent;
#endif
        // # proc per group
        extern int iNumProcPerGroup;
        
        
        //////// data type ////////
#ifndef _SERIAL_BUILD
        // datatype covertion, e.g., int -> MPI_INT
        template <typename T>
        inline MPI_Datatype typeMPI() = delete;
        
        // specializations
        template <>
        inline MPI_Datatype typeMPI<int>() {
            return MPI_INT;
        }
        
        template <>
        inline MPI_Datatype typeMPI<float>() {
            return MPI_FLOAT;
        }
        
        template <>
        inline MPI_Datatype typeMPI<double>() {
            return MPI_DOUBLE;
        }
        
        template <>
        inline MPI_Datatype typeMPI<std::complex<float>>() {
            return MPI_CXX_FLOAT_COMPLEX;
        }
        
        template <>
        inline MPI_Datatype typeMPI<std::complex<double>>() {
            return MPI_CXX_DOUBLE_COMPLEX;
        }
        
        template <>
        inline MPI_Datatype typeMPI<char>() {
            return MPI_CHAR;
        }
#endif
    }
    
    
    ////////////////////////////// basics //////////////////////////////
    // initialize MPI
    void initialize(int argc, char * argv[]);
    
    // finalize MPI
    void finalize();
    
    // setup group
    void setupGroup(int nprocPerGroup);
    
    // free group MPI_Comm
    void freeGroupComm();
    
    // number of processors
    inline int nproc() {
#ifndef _SERIAL_BUILD
        static int nproc;
        MPI_Comm_size(internal::iCommCurrent, &nproc);
        return nproc;
#else
        return 1;
#endif
    }
    
    // MPI rank
    inline int rank() {
#ifndef _SERIAL_BUILD
        static int rank;
        MPI_Comm_rank(internal::iCommCurrent, &rank);
        return rank;
#else
        return 0;
#endif
    }
    
    // MPI rank string
    inline std::string strRank() {
        std::stringstream ss;
        ss << rank();
        return ss.str();
    }
    
    // change current group level
    inline void enterWorld() {
#ifndef _SERIAL_BUILD
        internal::iCommCurrent = MPI_COMM_WORLD;
#endif
    }
    
    // enter super
    inline void enterSuper() {
#ifndef _SERIAL_BUILD
        internal::iCommCurrent = internal::iCommSuper;
#endif
    }
    
    // enter infer
    inline void enterInfer() {
#ifndef _SERIAL_BUILD
        internal::iCommCurrent = internal::iCommInfer;
#endif
    }
    
    // root
    inline bool root() {
        return rank() == 0;
    }
    
    // super
    // Am I a super proc in my group?
    // Must be called in World
    inline bool super() {
        return rank() % internal::iNumProcPerGroup == 0;
    }
    
    // number of groups
    // Must be called in World
    inline int ngroup() {
        return (nproc() - 1) / internal::iNumProcPerGroup + 1;
    }
    
    // barrier
    void barrier();
    
    // abort
    void abort(int err = 1);
    
    
    ////////////////////////////// broadcast //////////////////////////////
    // raw array, ALLOCATED
    template <typename T>
    void bcast(T *valptr, int size, int src = 0) {
#ifndef _SERIAL_BUILD
        MPI_Bcast(valptr, size, internal::typeMPI<T>(), src,
                  internal::iCommCurrent);
#endif
    }
    
    // single
    template <typename T>
    void bcast(T &value, int src = 0) {
        bcast(&value, 1, src);
    }
    
    // std::vector
    template <typename T>
    void bcast(std::vector<T> &vec, int src = 0) {
        // size
        int size = 0;
        if (rank() == src) {
            size = (int)vec.size();
        }
        bcast(size, src);
        // allocate
        if (rank() != src) {
            vec.resize(size);
        }
        // data
        bcast(vec.data(), size, src);
    }
    
    // Eigen::Matrix
    template <typename EigenMat>
    void bcastEigen(EigenMat &mat, int src = 0) {
        // size
        int dim[2];
        if (rank() == src) {
            dim[0] = (int)mat.rows();
            dim[1] = (int)mat.cols();
        }
        bcast(dim, 2, src);
        // allocate
        if (rank() != src) {
            mat.resize(dim[0], dim[1]);
        }
        // data
        bcast(mat.data(), (int)mat.size(), src);
    }
    
    // specialization for string
    void bcast(std::string &str, int src = 0);
    
    // specialization for vector of string
    void bcast(std::vector<std::string> &vecStr, int src = 0);
    
    
    ////////////////////////////// isend / irecv //////////////////////////////
    // isend, for ALLOCATED Eigen::Matrix
    template <typename EigenMat>
    void isend(int dest, const EigenMat &mat, MPI_Request &request) {
#ifndef _SERIAL_BUILD
        MPI_Isend(mat.data(), (int)mat.size(),
                  internal::typeMPI<typename EigenMat::Scalar>(),
                  dest, dest, internal::iCommCurrent, &request);
#endif
    }
    
    // irecv, for ALLOCATED Eigen::Matrix
    template <typename EigenMat>
    void irecv(int source, EigenMat &mat, MPI_Request &request) {
#ifndef _SERIAL_BUILD
        MPI_Irecv(mat.data(), (int)mat.size(),
                  internal::typeMPI<typename EigenMat::Scalar>(),
                  source, rank(), internal::iCommCurrent, &request);
#endif
    }
    
    // wait_all: must be implemented in .cpp
    void wait_all(int count, MPI_Request requests[]);
    
    
    ////////////////////////////// reduce //////////////////////////////
    // min
    template <typename T>
    T min(const T &value) {
        T result = value;
#ifndef _SERIAL_BUILD
        MPI_Allreduce(&value, &result, 1, internal::typeMPI<T>(), MPI_MIN,
                      internal::iCommCurrent);
#endif
        return result;
    }
    
    // max
    template <typename T>
    T max(const T &value) {
        T result = value;
#ifndef _SERIAL_BUILD
        MPI_Allreduce(&value, &result, 1, internal::typeMPI<T>(), MPI_MAX,
                      internal::iCommCurrent);
#endif
        return result;
    }
    
    // sum
    template <typename T>
    T sum(const T &value) {
        T result = value;
#ifndef _SERIAL_BUILD
        MPI_Allreduce(&value, &result, 1, internal::typeMPI<T>(), MPI_SUM,
                      internal::iCommCurrent);
#endif
        return result;
    }
    
    // sum eigen in-place
    template <typename EigenMat>
    void sumEigen(EigenMat &value) {
        EigenMat total(value);
#ifndef _SERIAL_BUILD
        MPI_Allreduce(value.data(), total.data(), (int)value.size(),
                      internal::typeMPI<typename EigenMat::Scalar>(), MPI_SUM,
                      internal::iCommCurrent);
#endif
        value = total;
    }
    
    
    ////////////////////////////// gather //////////////////////////////
    // gather single
    template <typename T>
    void gather(T val, std::vector<T> &vecVal, int dest) {
        if (dest < 0) {
            vecVal.resize(nproc());
#ifndef _SERIAL_BUILD
            MPI_Allgather(&val, 1, internal::typeMPI<T>(), vecVal.data(), 1,
                          internal::typeMPI<T>(), internal::iCommCurrent);
#else
            vecVal[0] = val;
#endif
        } else {
            if (rank() == dest) {
                vecVal.resize(nproc());
            }
#ifndef _SERIAL_BUILD
            MPI_Gather(&val, 1, internal::typeMPI<T>(), vecVal.data(), 1,
                       internal::typeMPI<T>(), dest, internal::iCommCurrent);
#else
            vecVal[0] = val;
#endif
        }
    }
    
    // vector to vector<vector>
    template <typename T>
    void gather(const std::vector<T> &vec,
                std::vector<std::vector<T>> &vecVec, int dest) {
        // size
        int size = (int)vec.size();
        std::vector<int> vSize;
        gather(size, vSize, dest);
        
        // allocate flattened
        std::vector<int> vDisp;
        std::vector<T> vecVecFlatten;
        if (dest < 0 || dest == rank()) {
            int totalSize = 0;
            vDisp = std::vector<int>(nproc(), 0);
            for (int iproc = 0; iproc < nproc(); iproc++) {
                vDisp[iproc] = totalSize;
                totalSize += vSize[iproc];
            }
            vecVecFlatten.resize(totalSize);
        }
        
        // bcast flattened
        if (dest < 0) {
#ifndef _SERIAL_BUILD
            MPI_Allgatherv(vec.data(), size, internal::typeMPI<T>(),
                           vecVecFlatten.data(), vSize.data(), vDisp.data(),
                           internal::typeMPI<T>(), internal::iCommCurrent);
#else
            vecVecFlatten = vec;
#endif
        } else {
#ifndef _SERIAL_BUILD
            MPI_Gatherv(vec.data(), size, internal::typeMPI<T>(),
                        vecVecFlatten.data(), vSize.data(), vDisp.data(),
                        internal::typeMPI<T>(), dest, internal::iCommCurrent);
#else
            vecVecFlatten = vec;
#endif
        }
        
        // cast back to nested vector
        if (dest < 0 || dest == rank()) {
            vecVec.clear();
            vecVec.reserve(nproc());
            for (int iproc = 0; iproc < nproc(); iproc++) {
                std::vector<T> vecRank(vecVecFlatten.begin() + vDisp[iproc],
                                       vecVecFlatten.begin() + vDisp[iproc] +
                                       vSize[iproc]);
                vecVec.push_back(vecRank);
            }
        }
    }
    
    // specialization for string
    void gather(const std::string &str,
                std::vector<std::string> &vecStr, int dest);
    
    // specialization for vector of string
    void gather(const std::vector<std::string> &vecStr,
                std::vector<std::vector<std::string>> &vecVecStr, int dest);
    
    
    ////////////////////////////// scatter //////////////////////////////
    // single
    template <typename T>
    void scatter(const std::vector<T> &vecVal, T &val, int src) {
#ifndef _SERIAL_BUILD
        MPI_Scatter(vecVal.data(), 1, internal::typeMPI<T>(), &val, 1,
                    internal::typeMPI<T>(), src, internal::iCommCurrent);
#else
        val = vecVal[0];
#endif
    }
    
    // vector<vector> to vector
    template <typename T>
    void scatter(const std::vector<std::vector<T>> &vecVec,
                 std::vector<T> &vec, int src) {
        // to be formed on src
        std::vector<int> vSize;
        std::vector<int> vDisp;
        std::vector<T> vecVecFlatten;
        if (rank() == src) {
            int totalSize = 0;
            vSize = std::vector<int>(nproc(), 0);
            vDisp = std::vector<int>(nproc(), 0);
            for (int iproc = 0; iproc < nproc();iproc++) {
                vSize[iproc] = (int)vecVec[iproc].size();
                vDisp[iproc] = totalSize;
                totalSize += vSize[iproc];
                vecVecFlatten.insert(vecVecFlatten.end(),
                                     vecVec[iproc].begin(),
                                     vecVec[iproc].end());
            }
        }
        
        // scatter size
        int size = 0;
        scatter(vSize, size, src);
        vec.resize(size);
        
        // scatter flattened
#ifndef _SERIAL_BUILD
        MPI_Scatterv(vecVecFlatten.data(), vSize.data(), vDisp.data(),
                     internal::typeMPI<T>(), vec.data(), size,
                     internal::typeMPI<T>(), src, internal::iCommCurrent);
#else
        vec = vecVec[0];
#endif
    }
    
    
    ////////////////////////////// map //////////////////////////////
    // gather map
    template<typename T>
    void gather(const std::map<std::string, T> &map,
                std::vector<std::map<std::string, T>> &vecMap, int dest) {
        // seperate keys and vals
        std::vector<std::string> keys;
        std::vector<T> vals;
        for (auto it = map.begin(); it != map.end(); it++) {
            keys.push_back(it->first);
            vals.push_back(it->second);
        }
        
        // gather keys and vals
        std::vector<std::vector<std::string>> allKeys;
        std::vector<std::vector<T>> allVals;
        gather(keys, allKeys, dest);
        gather(vals, allVals, dest);
        
        // back to map
        if (dest < 0 || dest == rank()) {
            vecMap.clear();
            for (int iproc = 0; iproc < nproc(); iproc++) {
                std::map<std::string, T> mapIP;
                for (int ikey = 0; ikey < allKeys[iproc].size(); ikey++) {
                    mapIP.insert({allKeys[iproc][ikey], allVals[iproc][ikey]});
                }
                vecMap.push_back(mapIP);
            }
        }
    }
    
    // aggregate map by sum
    template <typename T>
    void aggregate(std::map<std::string, T> &map, int dest) {
        // gather
        std::vector<std::map<std::string, T>> allMaps;
        gather(map, allMaps, dest);
        // aggregate all into one map
        if (dest < 0 || dest == rank()) {
            map.clear();
            for (int iproc = 0; iproc < nproc(); iproc++) {
                for (auto it = allMaps[iproc].begin();
                     it != allMaps[iproc].end(); it++) {
                    map.insert({it->first, (T)0});
                    map.at(it->first) += it->second;
                }
            }
        }
    }
    
    
    ////////////////////////////// external //////////////////////////////
    // verbose
    std::string verbose();
}

#endif /* mpi_hpp */
