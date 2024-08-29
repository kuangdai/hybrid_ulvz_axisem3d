//
//  sbstream.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/5/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  smart streams based on SimpleBinStream
//  https://github.com/shaovoon/simplebinstream

#ifndef sbstream_hpp
#define sbstream_hpp

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Weverything"
#include "simplebinstream/SimpleBinStream.h"
#pragma clang diagnostic pop

#include "eigen.hpp"
#include <type_traits>

namespace sbs {
    using namespace simple;
    
    ///////////////// ofstream /////////////////
    class ofstream {
    public:
        // constructor
        ofstream(const std::string &fname) {
            // check file using std::ofstream because sbstream does not throw
            std::ofstream ofs(fname, std::fstream::binary);
            if (!ofs) {
                throw std::runtime_error("sbs::ofstream::ofstream || "
                                         "Error opening or creating output "
                                         "file: || " + fname);
            }
            ofs.close();
            // open
            mOFS.open(fname.c_str());
        }
        
        // close
        void close() {
            mOFS.close();
        }
        
        // scalar (disabled for eigen)
        template<typename T>
        typename std::enable_if<!std::is_base_of<Eigen::DenseBase<T>, T>::value,
        ofstream>::type &operator<<(const T &val) {
            try {
                mOFS << val;
                return *this;
            } catch (...) {
                throw std::runtime_error("sbs::ofstream::operator<< || "
                                         "IO error from SBStream.");
            }
        }
        
        // eigen
        template<typename T>
        typename std::enable_if<std::is_base_of<Eigen::DenseBase<T>, T>::value,
        ofstream>::type &operator<<(const Eigen::DenseBase<T> &val) {
            try {
                mOFS << (int)val.rows() << (int)val.cols();
                if (val.size() == 0) {
                    return *this;
                }
                mOFS.write((const char *)(&val(0, 0)), val.size() *
                           sizeof(typename Eigen::DenseBase<T>::Scalar));
                return *this;
            } catch (...) {
                throw std::runtime_error("sbs::ofstream::operator<< || "
                                         "IO error from SBStream.");
            }
        }
        
        // vector of scalar (enabled only for fundamental)
        template<typename T>
        typename std::enable_if<std::is_fundamental<T>::value,
        ofstream>::type &operator<<(const std::vector<T> &val) {
            try {
                mOFS << (int)val.size();
                if (val.size() == 0) {
                    return *this;
                }
                mOFS.write((const char *)(&val[0]), val.size() * sizeof(T));
                return *this;
            } catch (...) {
                throw std::runtime_error("sbs::ofstream::operator<< || "
                                         "IO error from SBStream.");
            }
        }
        
        // vector of non-fundamental such as eigen and string
        template<typename T>
        typename std::enable_if<!std::is_fundamental<T>::value,
        ofstream>::type &operator<<(const std::vector<T> &val) {
            try {
                mOFS << (int)val.size();
                for (const T &mat: val) {
                    (*this) << mat;
                }
                return *this;
            } catch (...) {
                throw std::runtime_error("sbs::ofstream::operator<< || "
                                         "IO error from SBStream.");
            }
        }
        
    private:
        // stream
        file_ostream<std::true_type> mOFS;
    };
    
    ///////////////// ifstream /////////////////
    class ifstream {
    public:
        // constructor
        ifstream(const std::string &fname) {
            // check file using std::ifstream because sbstream does not throw
            std::ifstream ifs(fname, std::fstream::binary);
            if (!ifs) {
                throw std::runtime_error("sbs::ifstream::ifstream || "
                                         "Error opening input file: || "
                                         + fname);
            }
            ifs.close();
            // open
            mIFS.open(fname.c_str());
        }
        
        // close
        void close() {
            mIFS.close();
        }
        
        // scalar (disabled for eigen)
        template<typename T>
        typename std::enable_if<!std::is_base_of<Eigen::DenseBase<T>, T>::value,
        ifstream>::type &operator>>(T &val) {
            try {
                mIFS >> val;
                return *this;
            } catch (...) {
                throw std::runtime_error("sbs::ifstream::operator>> || "
                                         "IO error from SBStream.");
            }
        }
        
        // eigen
        template<typename T>
        typename std::enable_if<std::is_base_of<Eigen::DenseBase<T>, T>::value,
        ifstream>::type &operator>>(Eigen::DenseBase<T> &val) {
            try {
                int row, col;
                mIFS >> row >> col;
                val.derived().resize(row, col);
                if (val.size() == 0) {
                    return *this;
                }
                mIFS.read((char *)(&val(0, 0)), val.size() *
                          sizeof(typename Eigen::DenseBase<T>::Scalar));
                return *this;
            } catch (...) {
                throw std::runtime_error("sbs::ifstream::operator>> || "
                                         "IO error from SBStream.");
            }
        }
        
        // vector of scalar (enabled only for fundamental)
        template<typename T>
        typename std::enable_if<std::is_fundamental<T>::value,
        ifstream>::type &operator>>(std::vector<T> &val) {
            try {
                int size;
                mIFS >> size;
                val.resize(size);
                if (val.size() == 0) {
                    return *this;
                }
                mIFS.read((char *)(&val[0]), val.size() * sizeof(T));
                return *this;
            } catch (...) {
                throw std::runtime_error("sbs::ifstream::operator>> || "
                                         "IO error from SBStream.");
            }
        }
        
        // vector of non-fundamental such as eigen and string
        template<typename T>
        typename std::enable_if<!std::is_fundamental<T>::value,
        ifstream>::type &operator>>(std::vector<T> &val) {
            try {
                int size;
                mIFS >> size;
                val.resize(size);
                for (T &mat: val) {
                    (*this) >> mat;
                }
                return *this;
            } catch (...) {
                throw std::runtime_error("sbs::ifstream::operator>> || "
                                         "IO error from SBStream.");
            }
        }
        
        // get: slower than >> but easier to use
        template<typename T>
        T get() {
            try{
                T val;
                (*this) >> val;
                return val;
            } catch (...) {
                throw std::runtime_error("sbs::ifstream::get || "
                                         "IO error from SBStream.");
            }
        }
        
    private:
        // stream
        file_istream<std::true_type> mIFS;
    };
}

#endif /* sbstream_hpp */
