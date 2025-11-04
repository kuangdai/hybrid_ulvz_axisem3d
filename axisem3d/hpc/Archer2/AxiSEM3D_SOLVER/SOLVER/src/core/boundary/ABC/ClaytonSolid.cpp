//
//  ClaytonSolid.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/5/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  Clayton-Enquist ABC for solid points

#include "ClaytonSolid.hpp"

// stream
#include "Domain.hpp"
#include "ClaytonSolid1D.hpp"
#include "ClaytonSolid3D.hpp"

// parameter constructor
ClaytonSolid::ClaytonSolid(const std::shared_ptr<SolidPoint> &sp):
mSolidPoint(sp) {
    // nothing
}

// stream constructor
ClaytonSolid::ClaytonSolid(sbs::ifstream &ifs, const Domain &domain):
mSolidPoint(domain.getSolidPoint(ifs.get<int>())) {
    // nothing
}


////////////////////////////////////////
//////////////// static ////////////////
////////////////////////////////////////

// stream-based creator
std::unique_ptr<ClaytonSolid>
ClaytonSolid::createFromStream(sbs::ifstream &ifs, const Domain &domain) {
    // create pointer
    const std::string &clsName = ifs.get<std::string>();
    if (clsName == "ClaytonSolid1D") {
        return std::make_unique<ClaytonSolid1D>(ifs, domain);
    } else if (clsName == "ClaytonSolid3D") {
        return std::make_unique<ClaytonSolid3D>(ifs, domain);
    } else {
        throw std::runtime_error("ClaytonSolid::createFromStream || "
                                 "Unknown derived class of ClaytonSolid: "
                                 + clsName);
    }
}
