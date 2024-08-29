//
//  SolidFluidCoupling.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/30/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  solid-fluid boundary condition

#ifndef SolidFluidCoupling_hpp
#define SolidFluidCoupling_hpp

// point
#include "eigen_point.hpp"
#include <memory>
class SolidPoint;
class FluidPoint;

// stream
#include "sbstream.hpp"
class Domain;

class SolidFluidCoupling {
public:
    // parameter constructor
    SolidFluidCoupling(const std::shared_ptr<SolidPoint> &sp,
                       const std::shared_ptr<FluidPoint> &fp);
    
    // stream constructor
    SolidFluidCoupling(sbs::ifstream &ifs, const Domain &domain);
    
    // destructor
    virtual ~SolidFluidCoupling() = default;
    
    // get solid point
    const std::shared_ptr<SolidPoint> &getSolidPoint() const {
        return mSolidPoint;
    }
    
    // apply coupling
    void apply() const;
    
    
    ////////////////////////////// virtual //////////////////////////////
protected:
    // check compatibility
    virtual void checkCompatibility(int nr) const;
    
public:
    // solid => fluid
    virtual void coupleSolidToFluid(const eigen::CMatX3 &solidDispl,
                                    eigen::CColX &fluidStiff) const = 0;
    
    // fluid => solid
    virtual void coupleFluidToSolid(const eigen::CColX &fluidStiff,
                                    eigen::CMatX3 &solidStiff) const = 0;
    
protected:
    // coupled solid-fluid pair
    const std::shared_ptr<SolidPoint> mSolidPoint;
    const std::shared_ptr<FluidPoint> mFluidPoint;
    
    
    ////////////////////////////////////////
    //////////////// static ////////////////
    ////////////////////////////////////////
    
public:
    // stream-based creator
    static std::unique_ptr<SolidFluidCoupling>
    createFromStream(sbs::ifstream &ifs, const Domain &domain);
};

#endif /* SolidFluidCoupling_hpp */
