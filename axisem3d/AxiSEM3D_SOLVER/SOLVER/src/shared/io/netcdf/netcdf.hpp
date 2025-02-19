//
//  netcdf.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/21/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  netcdf header

#ifndef netcdf_hpp
#define netcdf_hpp

extern "C" {
#include <netcdf.h>
}

#ifdef _USE_PARALLEL_NETCDF
extern "C" {
#include <netcdf_par.h>
}
#endif

#include "eigen.hpp"
#include "eigen_tensor.hpp"
#include "bstring.hpp"
#include <type_traits>

namespace netcdf {
    ///////////////////////// datatype /////////////////////////
    // datatype covertion, e.g., int -> NC_INT
    template <typename T>
    inline nc_type typeNC() = delete;
    
    // specializations
    template <>
    inline nc_type typeNC<char>() {
        return NC_CHAR;
    }
    
    template <>
    inline nc_type typeNC<int>() {
        return NC_INT;
    }
    
    template <>
    inline nc_type typeNC<float>() {
        return NC_FLOAT;
    }
    
    template <>
    inline nc_type typeNC<double>() {
        return NC_DOUBLE;
    }
    
    
    ///////////////////////// error /////////////////////////
    // error handler
    inline void error(int retval, const std::string &funcName,
                      const std::string &fname) {
        if (retval != NC_NOERR) {
            throw std::
            runtime_error("netcdf::error || "
                          "Error in NetCDF function: " + funcName + " || "
                          "Error code = " + bstring::toString(retval) + " || "
                          "NetCDF file: " + fname);
        }
    }
    
    
    ///////////////////////// type and id /////////////////////////
    // get id
    inline int varID(int fid, const std::string &vname,
                     const std::string &fname) {
        int varid = -1;
        if (nc_inq_varid(fid, vname.c_str(), &varid) != NC_NOERR) {
            throw std::runtime_error("netcdf::getVariableID || "
                                     "Error finding variable. || "
                                     "Variable name: " + vname + " || "
                                     "NetCDF file: " + fname);
        }
        return varid;
    }
    
    // get id with datatype check
    inline int varTypeID(int fid, const std::string &vname,
                         const std::string &fname,
                         const nc_type &requiredTypeNC) {
        // get id
        int varid = varID(fid, vname, fname);
        // get nc datatype
        nc_type typeInFile;
        error(nc_inq_vartype(fid, varid, &typeInFile),
              "nc_inq_vartype", fname);
        // compare datatype
        if (typeInFile != requiredTypeNC) {
            throw std::runtime_error("netcdf::checkGetVariableID || "
                                     "Inconsistent C++ and netcdf datatype. || "
                                     "Variable name: " + vname + " || "
                                     "NetCDF file: " + fname);
        }
        return varid;
    }
    
    
    ///////////////////////// allowed containers /////////////////////////
    // vector
    template <class T>
    int getVariableID(int fid, const std::string &vname,
                      const std::string &fname,
                      const std::vector<T> &vec) {
        return varTypeID(fid, vname, fname, typeNC<T>());
    }
    
    // string
    inline int getVariableID(int fid, const std::string &vname,
                             const std::string &fname,
                             const std::string &str) {
        return varTypeID(fid, vname, fname, typeNC<char>());
    }
    
    // Eigen::Matrix, RowMajor or column vector
    template <class Derived>
    typename std::enable_if<Eigen::DenseBase<Derived>::IsRowMajor ||
    Eigen::DenseBase<Derived>::ColsAtCompileTime == 1, int>::type
    getVariableID(int fid, const std::string &vname,
                  const std::string &fname,
                  const Eigen::DenseBase<Derived> &mat) {
        return varTypeID(fid, vname, fname, typeNC<typename
                         Eigen::DenseBase<Derived>::Scalar>());
    }
    
    // Eigen::Tensor, RowMajor
    template <class Scalar, int R> int
    getVariableID(int fid, const std::string &vname,
                  const std::string &fname,
                  const Eigen::Tensor<Scalar, R, Eigen::RowMajor> &tensor) {
        return varTypeID(fid, vname, fname, typeNC<Scalar>());
    }
    
    // Eigen::Tensor, ColMajor of rank one
    template <class Scalar, int R>
    typename std::enable_if<R == 1, int>::type
    getVariableID(int fid, const std::string &vname,
                  const std::string &fname,
                  const Eigen::Tensor<Scalar, R, Eigen::ColMajor> &tensor) {
        return varTypeID(fid, vname, fname, typeNC<Scalar>());
    }
}

#endif /* netcdf_hpp */
