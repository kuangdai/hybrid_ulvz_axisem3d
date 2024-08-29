//
//  SolidFluidCoupling.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/31/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  solid-fluid boundary condition

#include "SolidFluidCoupling.hpp"

// point
#include "SolidPoint.hpp"
#include "FluidPoint.hpp"

// stream
#include "SolidFluidCoupling1D.hpp"
#include "SolidFluidCoupling3D.hpp"
#include "Domain.hpp"

// parameter constructor
SolidFluidCoupling::
SolidFluidCoupling(const std::shared_ptr<SolidPoint> &sp,
                   const std::shared_ptr<FluidPoint> &fp):
mSolidPoint(sp), mFluidPoint(fp) {
    if (mSolidPoint->getMeshTag() != mFluidPoint->getMeshTag()) {
        throw std::runtime_error("SolidFluidCoupling::SolidFluidCoupling || "
                                 "The coupled solid and fluid points have "
                                 "different mesh tags (positions).");
    }
}

// stream constructor
SolidFluidCoupling::
SolidFluidCoupling(sbs::ifstream &ifs, const Domain &domain):
mSolidPoint(domain.getSolidPoint(ifs.get<int>())),
mFluidPoint(domain.getFluidPoint(ifs.get<int>())) {
    if (mSolidPoint->getMeshTag() != mFluidPoint->getMeshTag()) {
        throw std::runtime_error("SolidFluidCoupling::SolidFluidCoupling || "
                                 "The coupled solid and fluid points have "
                                 "different mesh tags (positions).");
    }
}


////////////////////////////// virtual //////////////////////////////
// compute coupling
void SolidFluidCoupling::apply() const {
    // this order matters!
    coupleSolidToFluid(mSolidPoint->getFields().mDispl,
                       mFluidPoint->getFields().mStiff);
    coupleFluidToSolid(mFluidPoint->getFields().mStiff,
                       mSolidPoint->getFields().mStiff);
}


////////////////////////////// virtual //////////////////////////////
// check compatibility
void SolidFluidCoupling::checkCompatibility(int nr) const {
    if (mSolidPoint->getNr() != nr || mFluidPoint->getNr() != nr) {
        throw std::runtime_error("SolidFluidCoupling::checkCompatibility || "
                                 "Incompatible sizes.");
    }
}


////////////////////////////////////////
//////////////// static ////////////////
////////////////////////////////////////

// stream-based creator
std::unique_ptr<SolidFluidCoupling>SolidFluidCoupling::
createFromStream(sbs::ifstream &ifs, const Domain &domain) {
    // create pointer
    const std::string &clsName = ifs.get<std::string>();
    if (clsName == "SolidFluidCoupling1D") {
        return std::make_unique<SolidFluidCoupling1D>(ifs, domain);
    } else if (clsName == "SolidFluidCoupling3D") {
        return std::make_unique<SolidFluidCoupling3D>(ifs, domain);
    } else {
        throw std::runtime_error("SolidFluidCoupling::createFromStream || "
                                 "Unknown derived class of SolidFluidCoupling: "
                                 + clsName);
    }
}
