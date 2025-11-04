//
//  eigen_tools.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/11/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  eigen tools

#ifndef eigen_tools_hpp
#define eigen_tools_hpp

#include "eigen.hpp"
#include "bstring.hpp"

namespace eigen_tools {
    // memory info
    template <class EigenMat>
    std::string memoryInfo(const EigenMat &mat, const std::string &title,
                           double &addToMemGB, int count = 1) {
        // type name
        const std::string &tname =
        bstring::typeName<typename EigenMat::Scalar>();
        // scalar size
        size_t scalar = sizeof(typename EigenMat::Scalar);
        // total size
        double memGB = count * (mat.size() * scalar) / 1e9;
        // format
        std::stringstream ss;
        ss << title << ": ";
        if (count > 1) {
            ss << "matrix count = " << count << "; ";
        }
        ss << "dimensions = " << mat.rows() << " x " << mat.cols() << "; ";
        ss << "scalar type = " << tname << " (" << scalar << " bytes); ";
        ss << "memory = " << memGB << " GB";
        // add to
        addToMemGB += memGB;
        return ss.str();
    }
    
    // compute 2 * exp(I * alpha * phi)
    template <typename T,
    typename CColX = Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1>>
    void computeTwoExpIAlphaPhi(int nu_1, double phi, CColX &twoExpIAlphaPhi) {
        // no factor two on order zero
        twoExpIAlphaPhi(0) = 1.;
        // non-zero orders
        static const T two = 2.;
        std::complex<T> i_phi(0., phi);
        for (int alpha = 1; alpha < nu_1; alpha++) {
            twoExpIAlphaPhi(alpha) = two * exp((T)alpha * i_phi);
        }
    }
    
    // compute real value of a Fourier series at phi with precomputed exp
    template <typename T, int M,
    typename RRowM = Eigen::Matrix<T, 1, M>,
    typename CColX = Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1>,
    typename CMatXM = Eigen::Matrix<std::complex<T>, Eigen::Dynamic, M>>
    void computeFourierAtPhiExp(const CMatXM &coeffs, int nu_1,
                                const CColX &twoExpIAlphaPhi, RRowM &reals) {
        // sum over Fourier orders
        reals = (twoExpIAlphaPhi.topRows(nu_1).asDiagonal() *
                 coeffs.topRows(nu_1)).real().colwise().sum();
    }
    
    // compute real value of a Fourier series at phi
    template <typename T, int M,
    typename RRowM = Eigen::Matrix<T, 1, M>,
    typename CColX = Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1>,
    typename CMatXM = Eigen::Matrix<std::complex<T>, Eigen::Dynamic, M>>
    void computeFourierAtPhi(const CMatXM &coeffs, int nu_1, double phi,
                             RRowM &reals) {
        // 2 * exp(I * alpha * phi)
        CColX twoExpIAlphaPhi(nu_1);
        computeTwoExpIAlphaPhi<T>(nu_1, phi, twoExpIAlphaPhi);
        // sum over Fourier orders
        computeFourierAtPhiExp<T, M>(coeffs, nu_1, twoExpIAlphaPhi, reals);
    }
}

#endif /* eigen_tools_hpp */
