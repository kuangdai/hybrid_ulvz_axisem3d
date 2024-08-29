//
//  AxialBoundary.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/29/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  axial boundary condition

#include "AxialBoundary.hpp"

// point
#include "SolidPoint.hpp"
#include "FluidPoint.hpp"

// domain
#include "Domain.hpp"
#include "Messaging.hpp"

// add solid point
void AxialBoundary::addPoint(const std::shared_ptr<SolidPoint> &sp) {
    mSolidPoints.push_back(sp);
    
    // make sure mass is 1D
    if (sp->typeInfo().find("3D") != std::string::npos) {
        throw std::runtime_error("AxialBoundary::addPoint || "
                                 "An axial point must have a 1D mass");
    }
}

// add fluid point
void AxialBoundary::addPoint(const std::shared_ptr<FluidPoint> &fp) {
    mFluidPoints.push_back(fp);
    
    // make sure mass is 1D
    if (fp->typeInfo().find("3D") != std::string::npos) {
        throw std::runtime_error("AxialBoundary::addPoint || "
                                 "An axial point must have a 1D mass");
    }
}

// apply axial masking
void AxialBoundary::apply() const {
    static const numerical::ComplexR czero = 0.;
    static const numerical::ComplexR cJ = {0., 1.};
    static const numerical::Real half = .5;
    
    // solid
    for (const std::shared_ptr<SolidPoint> &sp: mSolidPoints) {
        eigen::CMatX3 &stiff = sp->getFields().mStiff;
        
        // alpha = 0
        stiff(0, 0) = czero;
        stiff(0, 1) = czero;
        
        // higher orders
        if (sp->getNu() >= 1) {
            // alpha = 1
            numerical::ComplexR s0 = stiff(1, 0);
            numerical::ComplexR s1 = stiff(1, 1);
            stiff(1, 0) = (s0 - cJ * s1) * half;
            stiff(1, 1) = (s1 + cJ * s0) * half;
            stiff(1, 2) = czero;
            // alpha > 1
            stiff.bottomRows(sp->getNu() - 1).setZero();
        }
    }
    
    // fluid
    for (const std::shared_ptr<FluidPoint> &fp: mFluidPoints) {
        fp->getFields().mStiff.bottomRows(fp->getNu()).setZero();
    }
}

// count info
void AxialBoundary::
countInfo(const Messaging &msg, int &solid, int &fluid) const {
    solid = 0;
    for (const auto &point: mSolidPoints) {
        if (!msg.pointInSmallerRank(point)) {
            solid++;
        }
    }
    fluid = 0;
    for (const auto &point: mFluidPoints) {
        if (!msg.pointInSmallerRank(point)) {
            fluid++;
        }
    }
}


////////////////////////////////////////
//////////////// static ////////////////
////////////////////////////////////////

// stream-based creator
std::unique_ptr<const AxialBoundary> AxialBoundary::
createFromStream(sbs::ifstream &ifs, const Domain &domain) {
    // collection
    std::unique_ptr<AxialBoundary> ab = std::make_unique<AxialBoundary>();
    
    // solid
    const std::vector<int> &spTags = ifs.get<std::vector<int>>();
    ab->mSolidPoints.reserve(spTags.size());
    for (int spTag: spTags) {
        ab->mSolidPoints.push_back(domain.getSolidPoint(spTag));
    }
    
    // fluid
    const std::vector<int> &fpTags = ifs.get<std::vector<int>>();
    ab->mFluidPoints.reserve(fpTags.size());
    for (int fpTag: fpTags) {
        ab->mFluidPoints.push_back(domain.getFluidPoint(fpTag));
    }
    return ab;
}
