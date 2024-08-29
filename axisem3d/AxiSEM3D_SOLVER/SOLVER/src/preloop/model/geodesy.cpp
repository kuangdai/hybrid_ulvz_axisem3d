//
//  geodesy.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/16/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  geodesy of the model and transformations between these systems:
//  Geographic: (lat, lon), reindeers at north pole
//  Geocentric: (theta, phi), reindeers at north pole
//  Source-centered: (distance, azimuth), source at north pole

#include "geodesy.hpp"
#include "bstring.hpp"
#include "InparamYAML.hpp"
#include "ExodusMesh.hpp"
#include "mpi.hpp"

namespace geodesy {
    // internal data
    namespace internal {
        // Cartesain
        bool iCartesian = false;
        // outer radius
        double iOuterRadius = -1;
        // depth zero
        double iDepthZeroRadius = -1.;
        // outer flattening
        double iOuterFlattening = -1.;
        // ellipticity
        std::vector<double> iEllipR;
        std::vector<double> iEllipF;
    }
    
    // setup
    void setup(const InparamYAML &inparamModel, ExodusMesh &exMesh) {
        // CS
        internal::iCartesian = exMesh.isCartesian();
        
        // radius
        double rMeshTop =
        inparamModel.getsWithOptions<double>("geodesy:radius_mesh_top", {
            {"AS_IS", exMesh.getMeshTop()}});
        // translation in z
        exMesh.translateZ(rMeshTop - exMesh.getMeshTop());
        // get outer radius after translation
        internal::iOuterRadius = exMesh.getMeshTop();
        
        // depth zero
        internal::iDepthZeroRadius =
        inparamModel.getsWithOptions<double>("geodesy:radius_depth_zero", {
            {"MESH_TOP", exMesh.getMeshTop()},
            {"SOLID_TOP", exMesh.getSolidTop()}});
        
        // flattening
        internal::iOuterFlattening =
        inparamModel.getsWithOptions<double>("geodesy:flattening", {
            {"SPHERE", 0.},
            {"WGS84", 1. / 298.257223563},
            {"GRS80", 1. / 298.257222101},
            {"SPECFEM3D_GLOBE", 1. / 299.8}});
        if (internal::iOuterFlattening < numerical::dEpsilon) {
            internal::iOuterFlattening = 0.;
        }
        
        // ellipticity (only on root in ExodusMesh)
        if (mpi::root()) {
            const eigen::DMatXX_RM &ellip = exMesh.getEllipticityCurve();
            internal::iEllipR.resize(ellip.cols());
            internal::iEllipF.resize(ellip.cols());
            for (int col = 0; col < ellip.cols(); col++) {
                internal::iEllipR[col] = ellip(0, col);
                internal::iEllipF[col] = ellip(1, col) *
                internal::iOuterFlattening;
            }
        }
        mpi::bcast(internal::iEllipR);
        mpi::bcast(internal::iEllipF);
    }
    
    // verbose
    std::string verbose(const eigen::DColX &discontinuities) {
        std::stringstream ss;
        ss << bstring::boxTitle("Geodesy");
        ss << bstring::boxEquals(0, 29, "model in Cartesian",
                                 internal::iCartesian);
        ss << bstring::boxEquals(0, 29, "outer radius",
                                 internal::iOuterRadius);
        ss << bstring::boxEquals(0, 29, "radius of depth zero",
                                 internal::iDepthZeroRadius);
        ss << bstring::boxEquals(0, 29, "flattening on surface",
                                 internal::iOuterFlattening);
        if (internal::iOuterFlattening > numerical::dEpsilon) {
            ///////// f at discontinuities /////////
            ss << bstring::boxSubTitle(0, "flattening at discontinuities");
            // compute f
            const eigen::DColX &fDisc = computeFlattening(discontinuities);
            // width
            int width = (int)bstring::toString(internal::iOuterFlattening *
                                               numerical::dPi / 3.14).size();
            // verbose f from surface to center
            for (int idisc = (int)fDisc.size() - 1; idisc >=0; idisc--) {
                ss << bstring::boxEquals(2, width + 6, "f = " +
                                         bstring::toString(fDisc[idisc]),
                                         discontinuities[idisc], "at");
            }
            ss << "Correcting for ellipticity in latitude-θ conversions.\n";
        } else {
            ss << "Perfect sphere, no ellipticity correction.\n";
        }
        ss << bstring::boxBaseline() << "\n\n";
        return ss.str();
    }
}
