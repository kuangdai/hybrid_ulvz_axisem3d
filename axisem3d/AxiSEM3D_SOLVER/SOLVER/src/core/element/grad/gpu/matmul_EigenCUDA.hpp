//
//  matmul_EigenCUDA.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 9/8/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  matrix multiplication on CUDA

#ifndef matmul_EigenCUDA_hpp
#define matmul_EigenCUDA_hpp

#include "eigen_element.hpp"

namespace matmul_EigenCUDA {
    // X = GT * F; Y = F * G
    void X_eq_GTxF__Y_eq_FxG(const eigen::RMatPP_RM &G,
                             const eigen::RMatPP_RM &GT,
                             const eigen::CMatPP_RM &F,
                             eigen::CMatPP_RM &X,
                             eigen::CMatPP_RM &Y);
    
    // F = G * X + Y * GT
    void F_eq_GxX_p_YxGT(const eigen::RMatPP_RM &G,
                         const eigen::RMatPP_RM &GT,
                         eigen::CMatPP_RM &F,
                         const eigen::CMatPP_RM &X,
                         const eigen::CMatPP_RM &Y);
}

#endif /* matmul_EigenCUDA_hpp */
