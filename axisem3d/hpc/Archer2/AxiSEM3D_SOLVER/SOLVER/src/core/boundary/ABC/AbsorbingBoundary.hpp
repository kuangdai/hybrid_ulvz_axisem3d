//
//  AbsorbingBoundary.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/1/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  absorbing boundary condition

#ifndef AbsorbingBoundary_hpp
#define AbsorbingBoundary_hpp

#include "ClaytonSolid.hpp"
#include "ClaytonFluid.hpp"
#include <vector>

// domain
#include <map>
class Messaging;

class AbsorbingBoundary {
public:
    // add Clayton Solid
    void addClaytonSolid(std::unique_ptr<const ClaytonSolid> &clayton) {
        mClaytonSolids.push_back(std::move(clayton));
    }
    
    // add Clayton Fluid
    void addClaytonFluid(std::unique_ptr<const ClaytonFluid> &clayton) {
        mClaytonFluids.push_back(std::move(clayton));
    }
    
    // apply ABC
    void apply() const {
        for (const std::unique_ptr<const ClaytonSolid> &clayton:
             mClaytonSolids) {
            clayton->apply();
        }
        for (const std::unique_ptr<const ClaytonFluid> &clayton:
             mClaytonFluids) {
            clayton->apply();
        }
    }
    
    // count info
    std::map<std::string, int> countInfo(const Messaging &msg) const;
    
private:
    // Claytons
    std::vector<std::unique_ptr<const ClaytonSolid>> mClaytonSolids;
    std::vector<std::unique_ptr<const ClaytonFluid>> mClaytonFluids;
    
    
    ////////////////////////////////////////
    //////////////// static ////////////////
    ////////////////////////////////////////
    
public:
    // stream-based creator
    static std::unique_ptr<const AbsorbingBoundary>
    createFromStream(sbs::ifstream &ifs, const Domain &domain);
};

#endif /* AbsorbingBoundary_hpp */
