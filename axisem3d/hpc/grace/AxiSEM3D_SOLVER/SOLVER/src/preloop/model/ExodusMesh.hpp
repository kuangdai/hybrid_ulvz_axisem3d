//
//  ExodusMesh.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 9/2/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  Exodus mesh created by salvus mesher

#ifndef ExodusMesh_hpp
#define ExodusMesh_hpp

#include "eigen_generic.hpp"
#include <map>

class InparamYAML;
class NetCDF_Reader;

namespace eigen {
    // connectivity
    typedef Eigen::Matrix<int, Eigen::Dynamic, 4, Eigen::RowMajor> IMatX4_RM;
    // coords
    typedef Eigen::Matrix<double, Eigen::Dynamic, 2> DMatX2;
}

class ExodusMesh {
public:
    // constructor
    ExodusMesh(const InparamYAML &inparamModel);
    
    // verbose
    std::string verbose() const;
    
    
    ///////////////// read and broadcast /////////////////
private:
    // global variables
    void readBcastGlobal(const NetCDF_Reader &reader,
                         double &memSup, double &memAll);
    // connectivity
    void readBcastConnectivity(const NetCDF_Reader &reader,
                               double &memSup, double &memAll);
    // coordinates
    void readBcastCoordinates(const NetCDF_Reader &reader,
                              double &memSup, double &memAll);
    // side sets
    void readBcastSideSets(const NetCDF_Reader &reader,
                           double &memSup, double &memAll);
    // elemental variables
    void readBcastElemental(const NetCDF_Reader &reader,
                            double &memSup, double &memAll);
    // radial variables
    void readBcastRadial(const NetCDF_Reader &reader,
                         double &memSup, double &memAll);
    // ellipticity
    void readBcastEllipticity(const NetCDF_Reader &reader,
                              double &memSup, double &memAll);
    
    
    ///////////////// data /////////////////
    // file name
    const std::string mFileName;
    
    // global variables and records
    std::map<std::string, double> mGlobalVariables;
    std::map<std::string, std::string> mGlobalRecords;
    // mesh generation cmdline
    std::string mCmdMeshGen;
    
    // connectivity (super only)
    eigen::IMatX4_RM mConnectivity;
    
    // coords (super only)
    eigen::DMatX2 mNodalCoords;
    
    // side sets
    std::string mKeyLeftSS;
    std::string mKeyRightSS;
    std::string mKeyBottomSS;
    std::string mKeyTopSS;
    std::map<std::string, std::map<int, int>> mSideSets;
    
    // element types (super only)
    eigen::IColX mElementTypes;
    
    // storage type
    bool mElementNodesStorage;
    
    // radial variables
    // NOTE: these are elemental variables depending ONLY on radius (depth)
    //       all material properties in Exodus are assumed to be radial
    eigen::DColX mRadialCoords;
    std::map<std::string, eigen::DColX> mRadialVariables;
    
    // discontinuities
    eigen::DColX mDiscontinuities;
    
    // ellipticity
    eigen::DMatXX_RM mEllipticityCurve;
    
    
public:
    ///////////////// coordinates /////////////////
    // mesh top
    double getMeshTop() const {
        return mNodalCoords.col(1).maxCoeff();
    }
    
    // solid top
    double getSolidTop() const {
        for (int irad = (int)mRadialCoords.size() - 1; irad > 0; irad--) {
            // get vs
            double vs = (isIsotropic() ?
                         mRadialVariables.at("VS")(irad) :
                         mRadialVariables.at("VSV")(irad));
            // solid
            if (vs > numerical::dEpsilon) {
                // shift top back by distTol
                return (mRadialCoords(irad) +
                        mGlobalVariables.at("min_edge_length") / 100.);
            }
        }
        throw std::runtime_error("ExodusMesh::getSolidTop || "
                                 "The model constains pure fluid.");
    }
    
    // vertical translation
    void translateZ(double dz) {
        mNodalCoords.array().col(1) += dz;
        mRadialCoords.array() += dz;
        mDiscontinuities.array() += dz;
        mEllipticityCurve.array().row(0) += dz;
    }
    

    ///////////////// get mesh properties /////////////////
    // number of nodes
    int getNumNodes() const {
        return (int)mNodalCoords.rows();
    }
    
    // number of quads
    int getNumQuads() const {
        return (int)mConnectivity.rows();
    }
    
    // Cartesian
    bool isCartesian() const {
        return mGlobalRecords.at("crdsys") != "spherical";
    }
    
    // attenuation
    bool hasAttenuation() const {
        return mGlobalVariables.find("nr_lin_solids") != mGlobalVariables.end();
    }
    
    // isotropic
    bool isIsotropic() const {
        return mRadialVariables.find("VP") != mRadialVariables.end();
    }
    
    // global
    double getGlobalVariable(const std::string &key) const {
        return mGlobalVariables.at(key);
    }
    
    // discontinuities
    const eigen::DColX &getDiscontinuities() const {
        return mDiscontinuities;
    }
    
    // ellipticity
    const eigen::DMatXX_RM &getEllipticityCurve() const {
        return mEllipticityCurve;
    }
    
    
    ///////////////// side methods /////////////////
    // get side
    int getSide(int iquad, const std::string &ssName) const {
        try {
            return mSideSets.at(ssName).at(iquad);
        } catch (...) {
            return -1;
        }
    }
    
    // left
    int getLeftSide(int iquad) const {
        return getSide(iquad, mKeyLeftSS);
    }
    
    // right
    int getRightSide(int iquad) const {
        return getSide(iquad, mKeyRightSS);
    }
    
    // bottom
    int getBottomSide(int iquad) const {
        return getSide(iquad, mKeyBottomSS);
    }
    
    // top
    int getTopSide(int iquad) const {
        return getSide(iquad, mKeyTopSS);
    }
    
    // axial
    int getAxialSide(int iquad) const {
        return getLeftSide(iquad);
    }
};

#endif /* ExodusMesh_hpp */
