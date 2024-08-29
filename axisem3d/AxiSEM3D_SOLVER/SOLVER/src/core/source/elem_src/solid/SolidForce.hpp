//
//  SolidForce.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/5/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  force source on solid element

#ifndef SolidForce_hpp
#define SolidForce_hpp

#include "SolidSource.hpp"
#include "numerical.hpp"
#include "SourceTimeFunction.hpp"

class SolidForce: public SolidSource {
public:
    // stf type
    typedef SourceTimeFunction<numerical::Real, 3> STF;
    
    // constructor
    SolidForce(const std::shared_ptr<const SolidElement> &element,
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
    inline static STF::CTMatXND sPattern = STF::CTMatXND(0, spectral::nPEM * 3);
};

#endif /* SolidForce_hpp */
