//
//  Element.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/1/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  spectral element

#include "Element.hpp"

// point
#include "Point.hpp"

// transform
#include "CoordTransformCartesian.hpp"
#include "CoordTransformSpherical.hpp"
#include "geodesy.hpp"

// stream
#include "Domain.hpp"
#include "SolidElement.hpp"
#include "FluidElement.hpp"

// side
#include "se_tools.hpp"

using spectral::nPEM;

/////////////////////////// point ///////////////////////////
// get point nu
std::array<int, spectral::nPEM> Element::getPointNu_1() const {
    std::array<int, spectral::nPEM> res;
    for (int ipnt = 0; ipnt < nPEM; ipnt++) {
        res[ipnt] = getPoint(ipnt).getNu() + 1;
    }
    return res;
}

// point set
void Element::pointSet(bool elemInFourier) {
    // order
    mNr = 0;
    for (int ipnt = 0; ipnt < nPEM; ipnt++) {
        mNr = std::max(getPoint(ipnt).getNr(), mNr);
    }
    mNu = mNr / 2;
    
    // check compatibility
    mGradQuad->checkCompatibility(mNr);
    if (mPRT) {
        mPRT->checkCompatibility(mNr, elemInFourier);
        // set to RTZ for PRT
        setToRTZ_ByMaterialOrPRT();
    }
}

// find boundary points by tag
std::vector<int> Element::
findBoundaryPointsByTag(const std::vector<int> &boundaryMeshTags) const {
    // point indices on an element boundary
    const std::vector<int> &allSidePoints = se_tools::getPointsOnAllSides();
    
    // result
    std::vector<int> pointsFound;
    pointsFound.reserve(allSidePoints.size());
    
    // check mesh tag
    for (int ipnt: allSidePoints) {
        if (vector_tools::foundSU(boundaryMeshTags,
                                  getPoint(ipnt).getMeshTag())) {
            // found a point on the injection boundary
            pointsFound.push_back(ipnt);
        }
    }
    
    // already sorted and unique
    return pointsFound;
}

// find boundary points by crds
std::vector<int> Element::
findBoundaryPointsByCrds(const std::vector<double> &boundaryCrdsRorZ,
                         const std::vector<double> &boundaryCrdsTorS,
                         double distTol) const {
    // point indices on an element boundary
    const std::vector<int> &allSidePoints = se_tools::getPointsOnAllSides();
    
    // collect coords
    eigen::DColX crdRorZ(nPEM);
    eigen::DColX crdTorS(nPEM);
    for (int ipnt = 0; ipnt < nPEM; ipnt++) {
        const eigen::DRow2 &sz = getPoint(ipnt).getCoords();
        if (geodesy::isCartesian()) {
            crdRorZ(ipnt) = sz(1);
            crdTorS(ipnt) = sz(0);
        } else {
            crdRorZ(ipnt) = sz.norm();
            crdTorS(ipnt) = acos(sz(1) / std::max(crdRorZ(ipnt),
                                                  numerical::dEpsilon));
        }
    }
    
    // result
    std::vector<int> pointsFound;
    
    // check r or z
    for (double bRorZ: boundaryCrdsRorZ) {
        for (int ipnt: allSidePoints) {
            if (std::abs(crdRorZ[ipnt] - bRorZ) < distTol) {
                pointsFound.push_back(ipnt);
            }
        }
    }
    
    // check t or s
    if (!geodesy::isCartesian()) {
        // change distance tolerance to angular tolerance
        distTol /= crdRorZ(nPEM / 2);
    }
    for (double bTorS: boundaryCrdsTorS) {
        for (int ipnt: allSidePoints) {
            if (std::abs(crdTorS[ipnt] - bTorS) < distTol) {
                pointsFound.push_back(ipnt);
            }
        }
    }
    
    // sort and unique points
    vector_tools::sortUnique(pointsFound);
    return pointsFound;
}

/////////////////////////// crd transform ///////////////////////////
// set to RTZ by material or PRT
void Element::setToRTZ_ByMaterialOrPRT() {
    // coordinate transform
    createCoordTransform();
}

// create coordinate transform
// internally called by setToRTZ_ByMaterialOrPRT
// externally called by prepareWavefieldOutput
void Element::createCoordTransform() {
    if (!mTransform) {
        // difference between Cartesian and spherical
        if (geodesy::isCartesian()) {
            mTransform = std::make_unique<const CoordTransformCartesian>();
        } else {
            // compute theta
            eigen::DMatPP_RM theta;
            for (int ipol = 0; ipol < spectral::nPED; ipol++) {
                for (int jpol = 0; jpol < spectral::nPED; jpol++) {
                    int ipnt = ipol * spectral::nPED + jpol;
                    const eigen::DRow2 &sz = getPoint(ipnt).getCoords();
                    double r = sz.norm();
                    theta(ipol, jpol) = r < numerical::dEpsilon ?
                    0. : acos(sz(1) / r);
                }
            }
            mTransform = std::make_unique<const CoordTransformSpherical>(theta);
        }
    }
}


////////////////////////////////////////
//////////////// static ////////////////
////////////////////////////////////////

// stream-based creator
std::shared_ptr<Element> Element::
createFromStream(sbs::ifstream &ifs, const Domain &domain) {
    const std::string &clsName = ifs.get<std::string>();
    if (clsName == "SolidElement") {
        return std::make_shared<SolidElement>(ifs, domain);
    } else if (clsName == "FluidElement") {
        return std::make_shared<FluidElement>(ifs, domain);
    } else {
        throw std::runtime_error("Element::createFromStream || "
                                 "Unknown derived class of Element: "
                                 + clsName);
    }
}
