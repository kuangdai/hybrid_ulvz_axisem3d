//
//  ClaytonFluid3D.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/30/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  Clayton-Enquist ABC for fluid points in 3D

#include "ClaytonFluid3D.hpp"
#include "FluidPoint.hpp"
#include "fftw.hpp"

// parameter constructor
ClaytonFluid3D::
ClaytonFluid3D(const std::shared_ptr<FluidPoint> &fp,
               const eigen::DColX &rho, const eigen::DColX &vp,
               const eigen::DColX &area):
ClaytonFluid(fp),
mAreaOverRhoVp(area.cwiseQuotient(rho.cwiseProduct(vp))
               .cast<numerical::Real>()) {
    // check compatibility
    checkCompatibility();
}

// stream constructor
ClaytonFluid3D::
ClaytonFluid3D(sbs::ifstream &ifs, const Domain &domain):
ClaytonFluid(ifs, domain),
mAreaOverRhoVp(ifs.get<eigen::RColX>()) {
    // check compatibility
    checkCompatibility();
}

// check compatibility
void ClaytonFluid3D::checkCompatibility() {
    // check size
    int nr = mFluidPoint->getNr();
    if (nr != mAreaOverRhoVp.rows()) {
        throw std::runtime_error("ClaytonFluid3D::checkCompatibility ||"
                                 "Incompatible sizes.");
    }
    
    // workspace
    if (sVecR.rows() < nr) {
        sVecR.resize(nr);
        sVecC.resize(nr / 2 + 1);
    }
    
    // report request to FFTW
    fftw::gFFT_1.addNR(nr);
}

// apply ABC
void ClaytonFluid3D::apply() const {
    // get fields
    const eigen::CColX &veloc = mFluidPoint->getFields().mVeloc;
    eigen::CColX &stiff = mFluidPoint->getFields().mStiff;
    
    // constants
    int nr = mFluidPoint->getNr();
    int nu_1 = nr / 2 + 1;
    
    // FFT: Fouier => cardinal
    fftw::gFFT_1.computeC2R(veloc, sVecR, nr);
    
    // multiply by area / (rho * vp) in cardinal space
    sVecR.topRows(nr).array() *= mAreaOverRhoVp.array();
    
    // FFT: cardinal => Fourier
    fftw::gFFT_1.computeR2C(sVecR, sVecC, nr);
    
    // subtract
    stiff -= sVecC.topRows(nu_1);
}
