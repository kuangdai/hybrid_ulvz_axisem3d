//
//  NetCDF_Reader.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/13/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  NetCDF reader

#ifndef NetCDF_Reader_hpp
#define NetCDF_Reader_hpp

#include "numerical.hpp"
#include "netcdf.hpp"
#include <vector>

class NetCDF_Reader {
public:
    // destructor
    ~NetCDF_Reader();
    
    
    ////////////////// file system //////////////////
    // open
    void open(const std::string &fname);
    
    // open parallel
    void openParallel(const std::string &fname);
    
    // close
    void close();
    
    // is open
    bool isOpen() const {
        return mFileName != "";
    }
    
    
    ///////////// workflow1: info -> allocation -> read /////////////
    // get variable id with type check
    template <class Container>
    int getVariableID(const std::string &vname, const Container &val) const {
        return netcdf::getVariableID(mFileID, vname, mFileName, val);
    }
    
    // get dimensions
    void getVariableDimensions(int varid,
                               std::vector<numerical::Int> &dims) const {
        // get number of dimensions
        int ndims = -1;
        netcdf::error(nc_inq_varndims(mFileID, varid, &ndims),
                      "nc_inq_varndims", mFileName);
        dims.resize(ndims);
        
        // get id of dimensions
        std::vector<int> id_dims(ndims, -1);
        netcdf::error(nc_inq_vardimid(mFileID, varid, id_dims.data()),
                      "nc_inq_vardimid", mFileName);
        
        // get dimensions
        for (int idim = 0; idim < ndims; idim++) {
            netcdf::error(nc_inq_dimlen(mFileID, id_dims[idim],
                                        (size_t *)(&dims[idim])),
                          "nc_inq_dimlen", mFileName);
        }
    }
    
    // read data from variable chunk (no id and datatype check)
    template <class Container>
    void readVariable(int varid, const std::string &vname, Container &val,
                      const std::vector<numerical::Int> &starts,
                      const std::vector<numerical::Int> &counts) const {
        // sizeof(numerical::Int) = sizeof(size_t), cast directly
        if (nc_get_vara(mFileID, varid, (size_t *)starts.data(),
                        (size_t *)counts.data(), val.data()) != NC_NOERR) {
            throw std::runtime_error("NetCDF_Reader::readVariable || "
                                     "Error reading variable from file. || "
                                     "Variable name: " + vname + " || "
                                     "NetCDF file: " + mFileName);
        }
    }
    
    
    ///////////// workflow2: read + allocation /////////////
    // only for Eigen::Matrix
    template <class Derived>
    void readEigen(const std::string &vname,
                   Eigen::DenseBase<Derived> &mat) const {
        // access variable with datatype check
        int varid = getVariableID(vname, mat);
        
        // get dimensions
        std::vector<numerical::Int> dims;
        getVariableDimensions(varid, dims);
        numerical::Int dimTotal =
        std::accumulate(dims.begin(), dims.end(), 1,
                        std::multiplies<numerical::Int>());
        
        // allocate matrix
        bool allocated = false;
        int rowsFix = Eigen::DenseBase<Derived>::RowsAtCompileTime;
        int colsFix = Eigen::DenseBase<Derived>::ColsAtCompileTime;
        if (rowsFix != Eigen::Dynamic && colsFix != Eigen::Dynamic) {
            // must be perfect match
            if (dims == std::vector<numerical::Int>({rowsFix, colsFix})) {
                allocated = true;
            }
        } else if (rowsFix != Eigen::Dynamic) {
            // keep row
            if (dimTotal % rowsFix == 0) {
                mat.derived().resize(rowsFix, dimTotal / rowsFix);
                allocated = true;
            }
        } else if (colsFix != Eigen::Dynamic) {
            // keep col
            if (dimTotal % colsFix == 0) {
                mat.derived().resize(dimTotal / colsFix, colsFix);
                allocated = true;
            }
        } else {
            // keep first
            mat.derived().resize(dims[0], dimTotal / dims[0]);
            allocated = true;
        }
        
        // get data
        if (allocated) {
            netcdf::error(nc_get_var(mFileID, varid, mat.derived().data()),
                          "nc_get_var", mFileName);
        } else {
            throw std::runtime_error("NetCDF_Reader::readEigen || "
                                     "Error interpreting variable dimensions.");
        }
    }
    
    // read string
    void readString(const std::string &vname,
                    std::vector<std::string> &dest) const;
    
    
    ////////////////// data //////////////////
private:
    // file ID
    int mFileID = -1;
    
    // file name
    std::string mFileName = "";
};

#endif /* NetCDF_Reader_hpp */
