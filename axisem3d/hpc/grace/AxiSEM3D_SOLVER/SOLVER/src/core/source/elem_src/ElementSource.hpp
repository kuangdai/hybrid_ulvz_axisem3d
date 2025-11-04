//
//  ElementSource.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/4/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  source on element

#ifndef ElementSource_hpp
#define ElementSource_hpp

class ElementSource {
public:
    // destructor
    virtual ~ElementSource() = default;
    
    // apply source at a time step
    virtual void apply(int timestep, double time) const = 0;
};

#endif /* ElementSource_hpp */
