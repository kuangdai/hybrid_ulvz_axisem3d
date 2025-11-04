//
//  FluidPressure.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/5/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  pressure source on fluid element

#include "FluidPressure.hpp"
#include "FluidElement.hpp"

// constructor
FluidPressure::FluidPressure(const std::shared_ptr<const FluidElement> &element,
                             std::unique_ptr<const STF> &stf):
FluidSource(element), mSTF(stf.release()) {
    // prepare for pressure source
    element->preparePressureSource();
    
    // workspace
    int nu_1 = mSTF->getPatternNu_1();
    if (sPattern.rows() < nu_1) {
        sPattern.resize(nu_1, spectral::nPEM);
    }
}

// apply source at a time step
void FluidPressure::apply(int timestep, double time) const {
    // compute pattern
    mSTF->getPatternAtTimeStep(timestep, time, sPattern);
    
    // apply to element
    mElement->addPressureSource(sPattern, mSTF->getPatternNu_1());
}
