//
//  ClaytonFluid.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/29/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  Clayton-Enquist ABC for fluid points

#ifndef ClaytonFluid_hpp
#define ClaytonFluid_hpp

// point
#include <memory>
class FluidPoint;

// stream
#include "sbstream.hpp"
class Domain;

class ClaytonFluid {
public:
    // parameter constructor
    ClaytonFluid(const std::shared_ptr<FluidPoint> &fp);
    
    // stream constructor
    ClaytonFluid(sbs::ifstream &ifs, const Domain &domain);
    
    // get point
    const std::shared_ptr<FluidPoint> &getPoint() const {
        return mFluidPoint;
    }
    
    // destructor
    virtual ~ClaytonFluid() = default;
    
    // apply ABC
    virtual void apply() const = 0;
    
protected:
    // point
    const std::shared_ptr<FluidPoint> mFluidPoint;
    

    ////////////////////////////////////////
    //////////////// static ////////////////
    ////////////////////////////////////////
    
public:
    // stream-based creator
    static std::unique_ptr<ClaytonFluid>
    createFromStream(sbs::ifstream &ifs, const Domain &domain);
};

#endif /* ClaytonFluid_hpp */
