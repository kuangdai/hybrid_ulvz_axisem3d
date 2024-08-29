//
//  SolidForce.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/5/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  force source on solid element

#include "SolidForce.hpp"
#include "SolidElement.hpp"

// constructor
SolidForce::SolidForce(const std::shared_ptr<const SolidElement> &element,
                       std::unique_ptr<const STF> &stf):
SolidSource(element), mSTF(stf.release()) {
    // prepare for force source
    element->prepareForceSource();
    
    // workspace
    int nu_1 = mSTF->getPatternNu_1();
    if (sPattern.rows() < nu_1) {
        sPattern.resize(nu_1, spectral::nPEM * 3);
    }
}

// apply source at a time step
void SolidForce::apply(int timestep, double time) const {
    // compute pattern
    mSTF->getPatternAtTimeStep(timestep, time, sPattern);
    
    // apply to element
    mElement->addForceSource(sPattern, mSTF->getPatternNu_1());
}
