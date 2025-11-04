//
//  NewmarkTimeScheme.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/30/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  Newmark time scheme

#ifndef NewmarkTimeScheme_hpp
#define NewmarkTimeScheme_hpp

#include "TimeScheme.hpp"

class NewmarkTimeScheme: public TimeScheme {
public:
    // parameter constructor
    NewmarkTimeScheme(int verboseInterval, int stabilityInterval):
    TimeScheme(verboseInterval, stabilityInterval) {
        // nothing
    }
    
    // stream constructor
    NewmarkTimeScheme(sbs::ifstream &ifs):
    TimeScheme(ifs) {
        // nothing
    }
    
    
    //////////////////// point ////////////////////
    // create fields on a solid point
    void createFields(SolidPoint &point) const;
    
    // create fields on a fluid point
    void createFields(FluidPoint &point) const;
    
    
    //////////////////// timeloop ////////////////////
public:
    // solve
    void solve() const;
    
private:
    // update fields on solid points
    void update(const std::vector<std::shared_ptr<SolidPoint>> &points) const;
    
    // update fields on fluid points
    void update(const std::vector<std::shared_ptr<FluidPoint>> &points) const;
};

#endif /* NewmarkTimeScheme_hpp */
