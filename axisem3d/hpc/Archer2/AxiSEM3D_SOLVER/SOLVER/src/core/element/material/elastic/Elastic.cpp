//
//  Elastic.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 4/28/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  elastic material

#include "Elastic.hpp"
#include "Isotropic.hpp"
#include "TransverselyIsotropic.hpp"
#include "Anisotropic.hpp"

////////////////////////////////////////
//////////////// static ////////////////
////////////////////////////////////////

// stream-based creator
std::unique_ptr<Elastic> Elastic::createFromStream(sbs::ifstream &ifs) {
    const std::string &clsName = ifs.get<std::string>();
    if (clsName == "Isotropic") {
        return std::make_unique<Isotropic>(ifs);
    } else if (clsName == "TransverselyIsotropic") {
        return std::make_unique<TransverselyIsotropic>(ifs);
    } else if (clsName == "Anisotropic") {
        return std::make_unique<Anisotropic>(ifs);
    } else {
        throw std::runtime_error("Elastic::createFromStream || "
                                 "Unknown derived class of Elastic: "
                                 + clsName);
    }
}
