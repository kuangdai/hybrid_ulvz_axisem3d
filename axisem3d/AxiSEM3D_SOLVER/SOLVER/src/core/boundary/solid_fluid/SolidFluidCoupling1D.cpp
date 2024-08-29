//
//  SolidFluidCoupling1D.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/30/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  solid-fluid boundary condition in 1D

#include "SolidFluidCoupling1D.hpp"
#include "SolidPoint.hpp"

// parameter constructor
SolidFluidCoupling1D::
SolidFluidCoupling1D(const std::shared_ptr<SolidPoint> &sp,
                     const std::shared_ptr<FluidPoint> &fp,
                     double ns_unassmb, double nz_unassmb,
                     double ns_assmb, double nz_assmb,
                     double massFluid):
SolidFluidCoupling(sp, fp),
mNormalS_UnassembledMPI(ns_unassmb),
mNormalZ_UnassembledMPI(nz_unassmb),
mNormalS_AssembledMPI_InvMassFluid(ns_assmb / massFluid),
mNormalZ_AssembledMPI_InvMassFluid(nz_assmb / massFluid) {
    checkCompatibility(mSolidPoint->getNr());
}

// stream constructor
SolidFluidCoupling1D::
SolidFluidCoupling1D(sbs::ifstream &ifs, const Domain &domain):
SolidFluidCoupling(ifs, domain),
mNormalS_UnassembledMPI(ifs.get<numerical::Real>()),
mNormalZ_UnassembledMPI(ifs.get<numerical::Real>()),
mNormalS_AssembledMPI_InvMassFluid(ifs.get<numerical::Real>()),
mNormalZ_AssembledMPI_InvMassFluid(ifs.get<numerical::Real>()) {
    checkCompatibility(mSolidPoint->getNr());
}

// solid => fluid
void SolidFluidCoupling1D::coupleSolidToFluid(const eigen::CMatX3 &solidDispl,
                                              eigen::CColX &fluidStiff) const {
    fluidStiff += mNormalS_UnassembledMPI * solidDispl.col(0);
    fluidStiff += mNormalZ_UnassembledMPI * solidDispl.col(2);
}

// fluid => solid
void SolidFluidCoupling1D::coupleFluidToSolid(const eigen::CColX &fluidStiff,
                                              eigen::CMatX3 &solidStiff) const {
    solidStiff.col(0) -= mNormalS_AssembledMPI_InvMassFluid * fluidStiff;
    solidStiff.col(2) -= mNormalZ_AssembledMPI_InvMassFluid * fluidStiff;
}
