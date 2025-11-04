//
//  MassOceanLoad1D.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/25/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  1D mass with ocean load
//  theta: angle between surface normal and z-axis

#include "MassOceanLoad1D.hpp"

// parameter constructor
MassOceanLoad1D::MassOceanLoad1D(double mass, double massOcean, double theta):
mIMH_CosT2_p_IMV_SinT2(1. / mass * cos(theta) * cos(theta) +
                       1. / (mass + massOcean) * sin(theta) * sin(theta)),
mIMH_SinT2_p_IMV_CosT2(1. / mass * sin(theta) * sin(theta) +
                       1. / (mass + massOcean) * cos(theta) * cos(theta)),
mIMV_m_IMH_x_CosT_SinT((1. / (mass + massOcean) - 1. / mass) *
                       cos(theta) * sin(theta)),
mIMH(1 / mass) {
    // nothing
}

// stream constructor
MassOceanLoad1D::MassOceanLoad1D(sbs::ifstream &ifs):
mIMH_CosT2_p_IMV_SinT2(ifs.get<numerical::Real>()),
mIMH_SinT2_p_IMV_CosT2(ifs.get<numerical::Real>()),
mIMV_m_IMH_x_CosT_SinT(ifs.get<numerical::Real>()),
mIMH(ifs.get<numerical::Real>()) {
    // nothing
}

// check compatibility
void MassOceanLoad1D::checkCompatibility(int nr, bool solid) const {
    // must on solid
    if (!solid) {
        throw std::runtime_error("MassOceanLoad1D::checkCompatibility || "
                                 "Incompatible types: "
                                 "ocean load on fluid point.");
    }
    
    // expand workspace if needed
    int nu_1 = nr / 2 + 1;
    if (sStiff3_col0.rows() < nu_1) {
        sStiff3_col0.resize(nu_1);
    }
}

// compute accel in-place for solid
void MassOceanLoad1D::computeAccel(eigen::CMatX3 &stiff3) const {
    // copy s
    int nu_1 = (int)stiff3.rows();
    sStiff3_col0.topRows(nu_1) = stiff3.col(0);
    
    // s, z
    stiff3.col(0) = (mIMH_CosT2_p_IMV_SinT2 * sStiff3_col0.topRows(nu_1) +
                     mIMV_m_IMH_x_CosT_SinT * stiff3.col(2));
    stiff3.col(2) = (mIMV_m_IMH_x_CosT_SinT * sStiff3_col0.topRows(nu_1) +
                     mIMH_SinT2_p_IMV_CosT2 * stiff3.col(2));
    
    // phi
    stiff3.col(1) *= mIMH;
}
