//
//  ClaytonSolid3D.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/29/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  Clayton-Enquist ABC for solid points in 3D
//  f = v.n n * rpa + (v - v.n n) * rsa
//    = rsa v + v.k k, with k = sqrt(rpa - rsa) n
//  f -> stifness force
//  v -> velocity
//  n -> unit normal of surface
//  rsa -> rho * vs * area
//  rpa -> rho * vp * area

#include "ClaytonSolid3D.hpp"
#include "SolidPoint.hpp"
#include "fftw.hpp"

// parameter constructor
ClaytonSolid3D::
ClaytonSolid3D(const std::shared_ptr<SolidPoint> &sp,
               const eigen::DColX &rho,
               const eigen::DColX &vs, const eigen::DColX &vp,
               const eigen::DColX &area,
               const eigen::DMatX3 &unitNormal):
ClaytonSolid(sp),
mRSA(rho.cwiseProduct(vs.cwiseProduct(area)).cast<numerical::Real>()),
mK((rho.cwiseProduct((vp - vs).cwiseProduct(area)).cwiseSqrt().asDiagonal()
    * unitNormal).cast<numerical::Real>()) {
    // check compatibility
    checkCompatibility();
}

// stream constructor
ClaytonSolid3D::
ClaytonSolid3D(sbs::ifstream &ifs, const Domain &domain):
ClaytonSolid(ifs, domain),
mRSA(ifs.get<eigen::RColX>()),
mK(ifs.get<eigen::RMatX3>()) {
    // check compatibility
    checkCompatibility();
}

// check compatibility
void ClaytonSolid3D::checkCompatibility() {
    // check size
    int nr = mSolidPoint->getNr();
    if (nr != mRSA.rows()) {
        throw std::runtime_error("ClaytonSolid3D::checkCompatibility ||"
                                 "Incompatible sizes.");
    }
    
    // workspace
    if (sVR.rows() < nr) {
        sVR.resize(nr, 3);
        sAR.resize(nr, 3);
        sAC.resize(nr / 2 + 1, 3);
    }
    
    // report request to FFTW
    fftw::gFFT_3.addNR(nr);
}

// apply ABC
void ClaytonSolid3D::apply() const {
    // get fields
    const eigen::CMatX3 &veloc = mSolidPoint->getFields().mVeloc;
    eigen::CMatX3 &stiff = mSolidPoint->getFields().mStiff;
    
    // constants
    int nr = (int)mRSA.rows();
    int nu_1 = nr / 2 + 1;
    
    // FFT: Fouier => cardinal
    fftw::gFFT_3.computeC2R(veloc, sVR, nr);
    
    // a = rsa V
    sAR.topRows(nr) = mRSA.asDiagonal() * sVR.topRows(nr);
    
    // a -= V.k k
    sAR.topRows(nr) -= (sVR.topRows(nr).cwiseProduct(mK)
                        .rowwise().sum().asDiagonal() * mK);
    
    // FFT: cardinal => Fouier
    fftw::gFFT_3.computeR2C(sAR, sAC, nr);
    
    // subtract
    stiff -= sAC.topRows(nu_1);
}
