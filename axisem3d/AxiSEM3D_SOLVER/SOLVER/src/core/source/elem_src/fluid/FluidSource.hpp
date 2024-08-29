//
//  FluidSource.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/4/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  source on fluid element

#ifndef FluidSource_hpp
#define FluidSource_hpp

#include "ElementSource.hpp"
#include <memory>
class FluidElement;

class FluidSource: public ElementSource {
public:
    // constructor
    FluidSource(const std::shared_ptr<const FluidElement> &element):
    mElement(element) {
        // nothing
    }
    
    // destructor
    virtual ~FluidSource() = default;
    
protected:
    // element pointer
    const std::shared_ptr<const FluidElement> mElement;
};

#endif /* FluidSource_hpp */
