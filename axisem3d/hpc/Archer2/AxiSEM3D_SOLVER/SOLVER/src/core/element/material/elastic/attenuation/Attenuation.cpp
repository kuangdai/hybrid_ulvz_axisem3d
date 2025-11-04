//
//  Attenuation.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 2/25/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  attenuation

#include "Attenuation.hpp"
#include "AttenuationFull.hpp"
#include "AttenuationCG4.hpp"

////////////////////////////////////////
//////////////// static ////////////////
////////////////////////////////////////

// stream-based creator
std::unique_ptr<Attenuation> Attenuation::createFromStream(sbs::ifstream &ifs) {
    const std::string &clsName = ifs.get<std::string>();
    if (clsName == "NULL") {
        return nullptr;
    } else if (clsName == "AttenuationFull") {
        return std::make_unique<AttenuationFull>(ifs);
    } else if (clsName == "AttenuationCG4") {
        return std::make_unique<AttenuationCG4>(ifs);
    } else {
        throw std::runtime_error("Attenuation::createFromStream || "
                                 "Unknown derived class of Attenuation: "
                                 + clsName);
    }
    return nullptr;
}
