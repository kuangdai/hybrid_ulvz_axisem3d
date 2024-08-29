//
//  FluidPressure.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/5/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  pressure source on fluid element

#ifndef FluidPressure_hpp
#define FluidPressure_hpp

#include "FluidSource.hpp"
#include "numerical.hpp"
#include "SourceTimeFunction.hpp"

class FluidPressure: public FluidSource {
public:
    // stf type
    typedef SourceTimeFunction<numerical::Real, 1> STF;
    
    // constructor
    FluidPressure(const std::shared_ptr<const FluidElement> &element,
                  std::unique_ptr<const STF> &stf);
    
    // apply source at a time step
    void apply(int timestep, double time) const;
    
private:
    // source pattern
    const std::unique_ptr<const STF> mSTF;
    
    
    ////////////////////////////////////////
    //////////////// static ////////////////
    ////////////////////////////////////////
    
    // workspace
    inline static STF::CTMatXND sPattern = STF::CTMatXND(0, spectral::nPEM);
};

#endif /* FluidPressure_hpp */
