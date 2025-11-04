//
//  FluidSurfaceBoundary.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/30/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  stress-free boundary condition on fluid

#include "FluidSurfaceBoundary.hpp"

// point
#include "FluidPoint.hpp"

// domain
#include "Domain.hpp"
#include "Messaging.hpp"

// apply stress-free boundary condition on fluid
void FluidSurfaceBoundary::apply() const {
    // pressure ≡ 0 or accel ≡ 0
    // so, veloc = disp = everything ≡ 0
    for (const std::shared_ptr<FluidPoint> &fp: mFluidPoints) {
        fp->getFields().mStiff.setZero();
    }
}

// count info
int FluidSurfaceBoundary::
countInfo(const Messaging &msg) const {
    int count = 0;
    for (const auto &point: mFluidPoints) {
        if (!msg.pointInSmallerRank(point)) {
            count++;
        }
    }
    return count;
}


////////////////////////////////////////
//////////////// static ////////////////
////////////////////////////////////////

// stream-based creator
std::unique_ptr<const FluidSurfaceBoundary> FluidSurfaceBoundary::
createFromStream(sbs::ifstream &ifs, const Domain &domain) {
    // collection
    std::unique_ptr<FluidSurfaceBoundary> fsb =
    std::make_unique<FluidSurfaceBoundary>();
    
    // points
    const std::vector<int> &fpTags = ifs.get<std::vector<int>>();
    fsb->mFluidPoints.reserve(fpTags.size());
    for (int fpTag: fpTags) {
        fsb->mFluidPoints.push_back(domain.getFluidPoint(fpTag));
    }
    return fsb;
}
