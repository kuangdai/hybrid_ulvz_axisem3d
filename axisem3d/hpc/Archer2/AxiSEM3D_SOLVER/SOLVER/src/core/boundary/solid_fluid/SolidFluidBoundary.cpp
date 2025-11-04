//
//  SolidFluidBoundary.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/5/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  solid-fluid boundary condition

#include "SolidFluidBoundary.hpp"

// domain
#include "Messaging.hpp"
#include "SolidPoint.hpp"
#include "bstring.hpp"
#include "vector_tools.hpp"

// count info
std::map<std::string, int> SolidFluidBoundary::
countInfo(const Messaging &msg) const {
    std::map<std::string, int> countMap;
    for (const auto &sfc: mSFCs) {
        if (!msg.pointInSmallerRank(sfc->getSolidPoint())) {
            vector_tools::aggregate(countMap, bstring::typeName(*sfc), 1);
        }
    }
    return countMap;
}


////////////////////////////////////////
//////////////// static ////////////////
////////////////////////////////////////

// stream-based creator
std::unique_ptr<const SolidFluidBoundary> SolidFluidBoundary::
createFromStream(sbs::ifstream &ifs, const Domain &domain) {
    // collection
    std::unique_ptr<SolidFluidBoundary> sfb =
    std::make_unique<SolidFluidBoundary>();
    
    // size
    int nsfc = ifs.get<int>();
    sfb->mSFCs.reserve(nsfc);
    
    // SFCs
    for (int isfc = 0; isfc < nsfc; isfc++) {
        sfb->mSFCs.push_back(SolidFluidCoupling::
                             createFromStream(ifs, domain));
    }
    return sfb;
}
