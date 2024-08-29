//
//  mapPPvsN.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/9/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  conversion between vec_arD_MatPP_RM and MatXND
//  where D = dimension of data

#ifndef mapPPvsN_hpp
#define mapPPvsN_hpp

#include "eigen.hpp"
#include "spectral.hpp"

namespace mapPPvsN {
    using spectral::nPEM;
    
    // PP -> N
    template <typename vec_arD_MatPP_RM, typename MatXND,
    typename RowN = Eigen::Matrix<typename MatXND::Scalar, 1, nPEM>>
    void PP2N(const vec_arD_MatPP_RM &vaMPP, MatXND &MatN, int nu_1) {
        int ndim = (int)(MatN.cols() / nPEM);
        for (int alpha = 0; alpha < nu_1; alpha++) {
            for (int idim = 0; idim < ndim; idim++) {
                MatN.block(alpha, idim * nPEM, 1, nPEM) =
                Eigen::Map<const RowN>(vaMPP[alpha][idim].data());
            }
        }
    }
    
    // N -> PP
    template <typename vec_arD_MatPP_RM, typename MatXND,
    typename RowN = Eigen::Matrix<typename MatXND::Scalar, 1, nPEM>>
    void N2PP(const MatXND &MatN, vec_arD_MatPP_RM &vaMPP, int nu_1) {
        int ndim = (int)(MatN.cols() / nPEM);
        for (int alpha = 0; alpha < nu_1; alpha++) {
            for (int idim = 0; idim < ndim; idim++) {
                Eigen::Map<RowN>(vaMPP[alpha][idim].data()) =
                MatN.block(alpha, idim * nPEM, 1, nPEM);
            }
        }
    }
}

#endif /* mapPPvsN_hpp */
