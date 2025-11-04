//
//  AxialBoundary.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/29/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  axial boundary condition

#ifndef AxialBoundary_hpp
#define AxialBoundary_hpp

// point
#include <vector>
#include <memory>
class SolidPoint;
class FluidPoint;

// domain
#include "sbstream.hpp"
class Domain;
class Messaging;

class AxialBoundary {
public:
    // add solid point
    void addPoint(const std::shared_ptr<SolidPoint> &sp);
    
    // add fluid point
    void addPoint(const std::shared_ptr<FluidPoint> &fp);
    
    // apply axial masking
    void apply() const;
    
    // count info
    void countInfo(const Messaging &msg, int &solid, int &fluid) const;
    
private:
    // points on axis
    std::vector<std::shared_ptr<SolidPoint>> mSolidPoints;
    std::vector<std::shared_ptr<FluidPoint>> mFluidPoints;
    
    
    ////////////////////////////////////////
    //////////////// static ////////////////
    ////////////////////////////////////////
    
public:
    // stream-based creator
    static std::unique_ptr<const AxialBoundary>
    createFromStream(sbs::ifstream &ifs, const Domain &domain);
};

#endif /* AxialBoundary_hpp */
