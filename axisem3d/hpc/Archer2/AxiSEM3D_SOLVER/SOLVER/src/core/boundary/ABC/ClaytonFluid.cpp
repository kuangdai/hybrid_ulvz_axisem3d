//
//  ClaytonFluid.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/5/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  Clayton-Enquist ABC for fluid points

#include "ClaytonFluid.hpp"

// stream
#include "Domain.hpp"
#include "ClaytonFluid1D.hpp"
#include "ClaytonFluid3D.hpp"

// parameter constructor
ClaytonFluid::ClaytonFluid(const std::shared_ptr<FluidPoint> &fp):
mFluidPoint(fp) {
    // nothing
}

// stream constructor
ClaytonFluid::ClaytonFluid(sbs::ifstream &ifs, const Domain &domain):
mFluidPoint(domain.getFluidPoint(ifs.get<int>())) {
    // nothing
}


////////////////////////////////////////
//////////////// static ////////////////
////////////////////////////////////////

// stream-based creator
std::unique_ptr<ClaytonFluid>
ClaytonFluid::createFromStream(sbs::ifstream &ifs, const Domain &domain) {
    // create pointer
    const std::string &clsName = ifs.get<std::string>();
    if (clsName == "ClaytonFluid1D") {
        return std::make_unique<ClaytonFluid1D>(ifs, domain);
    } else if (clsName == "ClaytonFluid3D") {
        return std::make_unique<ClaytonFluid3D>(ifs, domain);
    } else {
        throw std::runtime_error("ClaytonFluid::createFromStream || "
                                 "Unknown derived class of ClaytonFluid: "
                                 + clsName);
    }
}
