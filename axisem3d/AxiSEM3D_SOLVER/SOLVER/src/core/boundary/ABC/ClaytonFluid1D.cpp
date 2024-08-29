//
//  ClaytonFluid1D.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/30/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  Clayton-Enquist ABC for fluid points in 1D

#include "ClaytonFluid1D.hpp"
#include "FluidPoint.hpp"

// parameter constructor
ClaytonFluid1D::
ClaytonFluid1D(const std::shared_ptr<FluidPoint> &fp,
               double rho, double vp, double area):
ClaytonFluid(fp), mAreaOverRhoVp(area / (rho * vp)) {
    // nothing
}

// stream constructor
ClaytonFluid1D::
ClaytonFluid1D(sbs::ifstream &ifs, const Domain &domain):
ClaytonFluid(ifs, domain), mAreaOverRhoVp(ifs.get<numerical::Real>()) {
    // nothing
}

// apply ABC
void ClaytonFluid1D::apply() const {
    // get fields
    const eigen::CColX &veloc = mFluidPoint->getFields().mVeloc;
    eigen::CColX &stiff = mFluidPoint->getFields().mStiff;
    
    // apply
    stiff -= veloc * mAreaOverRhoVp;
}
