//
//  se_tools.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/23/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  spectral element tools

#ifndef se_tools_hpp
#define se_tools_hpp

#include "spectral.hpp"
#include "vector_tools.hpp"

namespace se_tools {
    // points on a side
    inline std::vector<int> getPointsOnSide(int side) {
        // mesh structure
        int ipol0 = 0, ipol1 = 0, jpol0 = 0, jpol1 = 0;
        if (side == 0) {
            ipol0 = 0;
            ipol1 = spectral::nPol;
            jpol0 = jpol1 = 0;
        } else if (side == 1) {
            ipol0 = ipol1 = spectral::nPol;
            jpol0 = 0;
            jpol1 = spectral::nPol;
        } else if (side == 2) {
            ipol0 = 0;
            ipol1 = spectral::nPol;
            jpol0 = jpol1 = spectral::nPol;
        } else if (side == 3) {
            ipol0 = ipol1 = 0;
            jpol0 = 0;
            jpol1 = spectral::nPol;
        } else {
            throw std::runtime_error("se_tools::getPointsOnSide || "
                                     "Side index must be in {0, 1, 2, 3}.");
        }
        // points on side
        std::vector<int> sidePnts(spectral::nPED, -1);
        int ipnt = 0;
        for (int ipol = ipol0; ipol <= ipol1; ipol++) {
            for (int jpol = jpol0; jpol <= jpol1; jpol++) {
                sidePnts[ipnt++] = ipol * spectral::nPED + jpol;
            }
        }
        return sidePnts;
    }
    
    // points on all sides
    inline const std::vector<int> &getPointsOnAllSides() {
        static std::vector<int> allSidePnts;
        // only form it once
        if (allSidePnts.size() == 0) {
            for (int side = 0; side < 4; side++) {
                const std::vector<int> &sps = getPointsOnSide(side);
                allSidePnts.insert(allSidePnts.end(), sps.begin(), sps.end());
            }
            vector_tools::sortUnique(allSidePnts);
        }
        return allSidePnts;
    }
}

#endif /* se_tools_hpp */
