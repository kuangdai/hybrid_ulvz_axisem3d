//
//  geodesy.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/16/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  geodesy of the model and transformations between these CSs:
//  * geographic: (lat, lon), north pole on top
//  * geocentric: (theta, phi), north pole on top
//  * source-centered: (distance, azimuth), source on top

#ifndef geodesy_hpp
#define geodesy_hpp

#include "eigen_generic.hpp"
#include "vector_tools.hpp"

class InparamYAML;
class ExodusMesh;

namespace geodesy {
    // internal data
    namespace internal {
        // Cartesain
        extern bool iCartesian;
        // outer radius
        extern double iOuterRadius;
        // depth zero
        extern double iDepthZeroRadius;
        // outer flattening
        extern double iOuterFlattening;
        // ellipticity
        extern std::vector<double> iEllipR;
        extern std::vector<double> iEllipF;
    }
    
    // get
    inline bool isCartesian() {
        return internal::iCartesian;
    }
    
    // get
    inline double getOuterRadius() {
        return internal::iOuterRadius;
    }
    
    // get
    inline double getDepthZeroRadius() {
        return internal::iDepthZeroRadius;
    }
    
    // get
    inline double getOuterFlattening() {
        return internal::iOuterFlattening;
    }
    
    // compute flattenings at given radii
    template <class Container>
    Container computeFlattening(const Container &r) {
        Container f = r;
        int index0, index1;
        double factor0, factor1;
        for (int ir = 0; ir < r.size(); ir++) {
            vector_tools::linearInterpSorted(internal::iEllipR, r[ir],
                                             index0, index1, factor0, factor1);
            f[ir] = (internal::iEllipF[index0] * factor0 +
                     internal::iEllipF[index1] * factor1);
        }
        return f;
    }
    
    // setup
    void setup(const InparamYAML &inparamModel, ExodusMesh &exMesh);
    
    // verbose
    std::string verbose(const eigen::DColX &discontinuities);
}

#endif /* geodesy_hpp */
