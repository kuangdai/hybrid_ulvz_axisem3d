//
//  Mass.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/5/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  mass for both solid and fluid

#include "Mass.hpp"
#include "Mass1D.hpp"
#include "Mass3D.hpp"
#include "MassOceanLoad1D.hpp"
#include "MassOceanLoad3D.hpp"

////////////////////////////////////////
//////////////// static ////////////////
////////////////////////////////////////

// stream-based creator
std::unique_ptr<Mass> Mass::createFromStream(sbs::ifstream &ifs) {
    const std::string &clsName = ifs.get<std::string>();
    if (clsName == "Mass1D") {
        return std::make_unique<Mass1D>(ifs);
    } else if (clsName == "Mass3D") {
        return std::make_unique<Mass3D>(ifs);
    } else if(clsName == "MassOceanLoad1D") {
        return std::make_unique<MassOceanLoad1D>(ifs);
    } else if(clsName == "MassOceanLoad3D") {
        return std::make_unique<MassOceanLoad3D>(ifs);
    } else {
        throw std::runtime_error("Mass::createFromStream || "
                                 "Unknown derived class of Mass: " + clsName);
    }
}
