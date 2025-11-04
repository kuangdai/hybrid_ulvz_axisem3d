//
//  InparamYAML.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/31/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  YAML parser for input parameters
//  based on mini-yaml:
//  https://github.com/jimmiebergmann/mini-yaml

#include "InparamYAML.hpp"
#include "mpi.hpp"

// parse
void InparamYAML::parse(const std::string &fname) {
    // read input file
    std::string contents;
    if (mpi::root()) {
        contents = bstring::readAll(fname, "InparamYAML::parse");
    }
    mpi::bcast(contents);
    
    // parse input
    try {
        Yaml::Parse(mRoot, contents);
    } catch (...) {
        throw std::runtime_error("InparamYAML::parse || "
                                 "Error parsing YAML input file: || " + fname);
    }
}

// verbose
std::string InparamYAML::verbose() const {
    // get contents
    std::string contents;
    Yaml::Serialize(mRoot, contents);
    
    // add indent
    contents = "  " + bstring::replace(contents, "\n", "\n  ");
    contents = contents.substr(0, contents.size() - 2);
    
    // add subtitle
    std::stringstream ss;
    ss << bstring::boxSubTitle(0, mName);
    ss << contents;
    return ss.str();
}
