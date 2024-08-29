//
//  ClaytonSolid1D.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/29/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  Clayton-Enquist ABC for solid points in 1D
//  theta: angle between surface normal and z-axis

#ifndef ClaytonSolid1D_hpp
#define ClaytonSolid1D_hpp

#include "ClaytonSolid.hpp"
#include "numerical.hpp"

class ClaytonSolid1D: public ClaytonSolid {
public:
    // parameter constructor
    ClaytonSolid1D(const std::shared_ptr<SolidPoint> &sp,
                   double rho, double vs, double vp, double area, double theta);
    
    // stream constructor
    ClaytonSolid1D(sbs::ifstream &ifs, const Domain &domain);
    
    // apply ABC
    void apply() const;
    
private:
    // RSA = rho * vs * area
    // RPA = rho * vp * area
    // RSA Cos[t]^2 + RPA Sin[t]^2
    const numerical::Real mRSA_CosT2_p_RPA_SinT2;
    // RSA Sin[t]^2 + RPA Cos[t]^2
    const numerical::Real mRSA_SinT2_p_RPA_CosT2;
    // (RPA - RSA) Cos[t] Sin[t]
    const numerical::Real mRPA_m_RSA_x_CosT_SinT;
    // RSA (for transverse component)
    const numerical::Real mRSA;
};

#endif /* ClaytonSolid1D_hpp */
