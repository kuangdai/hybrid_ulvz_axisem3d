//
//  AbsorbingBoundary.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/5/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  absorbing boundary condition

#include "AbsorbingBoundary.hpp"
#include "ClaytonSolid.hpp"
#include "ClaytonFluid.hpp"

// domain
#include "Messaging.hpp"
#include "SolidPoint.hpp"
#include "FluidPoint.hpp"
#include "bstring.hpp"
#include "vector_tools.hpp"

// count info
std::map<std::string, int> AbsorbingBoundary::
countInfo(const Messaging &msg) const {
    std::map<std::string, int> countMap;
    // solid
    for (const auto &clayton: mClaytonSolids) {
        if (!msg.pointInSmallerRank(clayton->getPoint())) {
            vector_tools::aggregate(countMap, bstring::typeName(*clayton), 1);
        }
    }
    // fluid
    for (const auto &clayton: mClaytonFluids) {
        if (!msg.pointInSmallerRank(clayton->getPoint())) {
            vector_tools::aggregate(countMap, bstring::typeName(*clayton), 1);
        }
    }
    return countMap;
}


////////////////////////////////////////
//////////////// static ////////////////
////////////////////////////////////////

// stream-based creator
std::unique_ptr<const AbsorbingBoundary> AbsorbingBoundary::
createFromStream(sbs::ifstream &ifs, const Domain &domain) {
    // collection
    std::unique_ptr<AbsorbingBoundary> ab =
    std::make_unique<AbsorbingBoundary>();
    
    // solid
    int sizeSolid = ifs.get<int>();
    ab->mClaytonSolids.reserve(sizeSolid);
    for (int iab = 0; iab < sizeSolid; iab++) {
        ab->mClaytonSolids.push_back(ClaytonSolid::
                                     createFromStream(ifs, domain));
    }
    
    // fluid
    int sizeFluid = ifs.get<int>();
    ab->mClaytonFluids.reserve(sizeFluid);
    for (int iab = 0; iab < sizeFluid; iab++) {
        ab->mClaytonFluids.push_back(ClaytonFluid::
                                     createFromStream(ifs, domain));
    }
    return ab;
}
