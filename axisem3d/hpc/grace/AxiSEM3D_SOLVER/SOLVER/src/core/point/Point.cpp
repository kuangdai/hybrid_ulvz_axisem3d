//
//  Point.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/24/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  GLL point

#include "Point.hpp"
#include "bstring.hpp"
#include "SolidPoint.hpp"
#include "FluidPoint.hpp"

// type info
std::string Point::typeInfo() const {
    return bstring::typeName(*this) + "$" + bstring::typeName(*mMass);
}


////////////////////////////////////////
//////////////// static ////////////////
////////////////////////////////////////

// stream-based creator
std::shared_ptr<Point> Point::
createFromStream(sbs::ifstream &ifs, const TimeScheme &timeScheme) {
    const std::string &clsName = ifs.get<std::string>();
    if (clsName == "SolidPoint") {
        return std::make_shared<SolidPoint>(ifs, timeScheme);
    } else if (clsName == "FluidPoint") {
        return std::make_shared<FluidPoint>(ifs, timeScheme);
    } else {
        throw std::runtime_error("Point::createFromStream || "
                                 "Unknown derived class of Point: " + clsName);
    }
}
