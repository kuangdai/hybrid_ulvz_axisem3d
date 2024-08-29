//
//  ClaytonSolid1D.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/29/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  Clayton-Enquist ABC for solid points in 1D
//  theta: angle between surface normal and z-axis

#include "ClaytonSolid1D.hpp"
#include "SolidPoint.hpp"

// parameter constructor
ClaytonSolid1D::
ClaytonSolid1D(const std::shared_ptr<SolidPoint> &sp,
               double rho, double vs, double vp, double area, double theta):
ClaytonSolid(sp),
mRSA_CosT2_p_RPA_SinT2(rho * vs * area * cos(theta) * cos(theta) +
                       rho * vp * area * sin(theta) * sin(theta)),
mRSA_SinT2_p_RPA_CosT2(rho * vs * area * sin(theta) * sin(theta) +
                       rho * vp * area * cos(theta) * cos(theta)),
mRPA_m_RSA_x_CosT_SinT((rho * vp * area - rho * vs * area) *
                       cos(theta) * sin(theta)),
mRSA(rho * vs * area) {
    // nothing
}

// stream constructor
ClaytonSolid1D::
ClaytonSolid1D(sbs::ifstream &ifs, const Domain &domain):
ClaytonSolid(ifs, domain),
mRSA_CosT2_p_RPA_SinT2(ifs.get<numerical::Real>()),
mRSA_SinT2_p_RPA_CosT2(ifs.get<numerical::Real>()),
mRPA_m_RSA_x_CosT_SinT(ifs.get<numerical::Real>()),
mRSA(ifs.get<numerical::Real>()) {
    // nothing
}

// apply ABC
void ClaytonSolid1D::apply() const {
    // get fields
    const eigen::CMatX3 &veloc = mSolidPoint->getFields().mVeloc;
    eigen::CMatX3 &stiff = mSolidPoint->getFields().mStiff;
    
    // s, z
    stiff.col(0) -= (mRSA_CosT2_p_RPA_SinT2 * veloc.col(0) +
                     mRPA_m_RSA_x_CosT_SinT * veloc.col(2));
    stiff.col(2) -= (mRPA_m_RSA_x_CosT_SinT * veloc.col(0) +
                     mRSA_SinT2_p_RPA_CosT2 * veloc.col(2));
    
    // phi
    stiff.col(1) -= mRSA * veloc.col(1);
}
