// Volumetric3D_EMC_ULVZ.h
// created by Kuangdai on 16-May-2017 
// genetral Volumetric3D model with IRIS-EMC format

#pragma once

#include "Volumetric3D.h"
#include "eigenp.h"

class Volumetric3D_EMC_ULVZ: public Volumetric3D {
public:

    void initialize();
    void initialize(const std::vector<std::string> &params);
    bool get3dProperties(double r, double theta, double phi, double rElemCenter,
        std::vector<MaterialProperty> &properties, 
        std::vector<MaterialRefType> &refTypes,
        std::vector<double> &values) const;
    bool get3dProperties(double r, double theta, double phi, double rElemCenter,
        double thetaElemCenter,
        std::vector<MaterialProperty> &properties, 
        std::vector<MaterialRefType> &refTypes,
        std::vector<double> &values) const;    
        
    std::string verbose() const;

    void setSourceLocation(double srcLat, double srcLon, double srcDep) {
        mSrcLat = srcLat;
        mSrcLon = srcLon;
        mSrcDep = srcDep;
    };
    
private:
    
    // file
    std::string mFileName;
    std::string mVarName;
    
    // property
    MaterialProperty mMaterialProp;
    MaterialRefType mReferenceType;
    
    // factor
    double mFactor = 1.0;
    
    // source-centered
    double mSrcLat = 0.;
    double mSrcLon = 0.;
    double mSrcDep = 0.;

    // consider vertical disc or not
    bool mVerticalDiscontinuities = true;
    
    // data
    std::vector<RDMatXX> mGridData;
    RDColX mGridDep;
    RDColX mGridLat;
    RDColX mGridLon;
    
    // special model flag
    // abs -- use absolute value of the perturbations
    // pow -- use power of the perturbations
    // to use these special model flags, reference type cannot be Absolute
    // for "abs" flag, the following factor is used to change the sign
    // for "pow" flag, the following factor specifies the power, e.g.,
    // 2.0 means "squared", which makes the model sharper  
    // 0.5 means "sqrt", which makes the model smoother
    // the pow flag keeps the sign and global absolute maximum of the perturbations
    std::string mModelFlag = "none";
    double mModelFlagFactor = 1.0;
};

