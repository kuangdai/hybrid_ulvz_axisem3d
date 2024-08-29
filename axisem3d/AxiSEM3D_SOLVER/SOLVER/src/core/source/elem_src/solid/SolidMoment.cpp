//
//  SolidMoment.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/5/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  moment source on solid element

#include "SolidMoment.hpp"
#include "SolidElement.hpp"

// constructor
SolidMoment::SolidMoment(const std::shared_ptr<const SolidElement> &element,
                         std::unique_ptr<const STF> &stf):
SolidSource(element), mSTF(stf.release()) {
    // prepare for moment source
    element->prepareMomentSource();

    // workspace
    int nu_1 = mSTF->getPatternNu_1();
    if (sPattern.rows() < nu_1) {
        sPattern.resize(nu_1, spectral::nPEM * 6);
    }
}

// apply source at a time step
void SolidMoment::apply(int timestep, double time) const {
    // compute pattern
    mSTF->getPatternAtTimeStep(timestep, time, sPattern);
    
    // apply to element
    mElement->addMomentSource(sPattern, mSTF->getPatternNu_1());
}
