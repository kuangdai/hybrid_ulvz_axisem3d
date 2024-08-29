//
//  SolidFluidBoundary.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/1/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  solid-fluid boundary condition

#ifndef SolidFluidBoundary_hpp
#define SolidFluidBoundary_hpp

#include "SolidFluidCoupling.hpp"
#include <vector>

// domain
#include <map>
class Messaging;

class SolidFluidBoundary {
public:
    // add solid-fluid coupling
    void addSolidFluidCoupling(std::unique_ptr<const SolidFluidCoupling> &sf) {
        mSFCs.push_back(std::move(sf));
    }
    
    // apply solid-fluid coupling
    void apply() const {
        for (const std::unique_ptr<const SolidFluidCoupling> &sfc: mSFCs) {
            sfc->apply();
        }
    }
    
    // count info
    std::map<std::string, int> countInfo(const Messaging &msg) const;
    
private:
    // solid-fluid coupling pairs
    std::vector<std::unique_ptr<const SolidFluidCoupling>> mSFCs;
    
    
    ////////////////////////////////////////
    //////////////// static ////////////////
    ////////////////////////////////////////
    
public:
    // stream-based creator
    static std::unique_ptr<const SolidFluidBoundary>
    createFromStream(sbs::ifstream &ifs, const Domain &domain);
};

#endif /* SolidFluidBoundary_hpp */
