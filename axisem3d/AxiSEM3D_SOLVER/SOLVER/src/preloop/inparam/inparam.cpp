//
//  inparam.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 9/1/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  global input parameters

#include "inparam.hpp"
#include "io.hpp"
#include "bstring.hpp"

namespace inparam {
    // global input parameters
    InparamYAML gInparamModel("Model");
    InparamYAML gInparamNu("Nu");
    InparamYAML gInparamTime("Time");
    InparamYAML gInparamSource("Source");
    InparamYAML gInparamOutput("Output");
    InparamYAML gInparamAdvanced("Advanced");
    
    // setup
    void setup() {
        // parse
        gInparamModel.parse(io::gInputDirectory + "/inparam.model.yaml");
        gInparamAdvanced.parse(io::gInputDirectory + "/inparam.advanced.yaml");
        
        // io verbose
        io::gVerbose = gInparamAdvanced.
        getsWithLimits<io::VerboseLevel>("verbose:level", {
            {"none", io::VerboseLevel::none},
            {"essential", io::VerboseLevel::essential},
            {"detailed", io::VerboseLevel::detailed}});
        io::gVerboseWarnings =
        gInparamAdvanced.gets<bool>("verbose:warnings");
    }
    
    // verbose
    std::string verbose() {
        std::stringstream ss;
        ss << bstring::boxTitle("Parameters");
        ss << gInparamModel.verbose();
        ss << gInparamNu.verbose();
        ss << gInparamTime.verbose();
        ss << gInparamSource.verbose();
        ss << gInparamOutput.verbose();
        ss << gInparamAdvanced.verbose();
        ss << bstring::boxBaseline() << "\n\n";
        return ss.str();
    }
}
