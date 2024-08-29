//
//  Element.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/1/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  spectral element

#ifndef Element_hpp
#define Element_hpp

// components
#include "CoordTransform.hpp"
#include "GradientQuadrature.hpp"
#include "PRT.hpp"

#include <array>

// Real-typed GradientQuadrature
typedef GradientQuadrature<numerical::Real> GradQuad;

// point
class Point;
class Domain;

class Element {
public:
    // parameter constructor
    Element(int quadTag, std::unique_ptr<const GradQuad> &gq,
            std::unique_ptr<const PRT> &prt):
    mQuadTag(quadTag), mGradQuad(gq.release()), mPRT(prt.release()) {
        mTransform = nullptr;
    }
    
    // copy constructor
    Element(const Element &other):
    mQuadTag(other.mQuadTag),
    mGradQuad(std::make_unique<GradQuad>(*(other.mGradQuad))),
    mPRT(other.mPRT == nullptr ? nullptr :
         std::make_unique<PRT>(*(other.mPRT))) {
        mTransform = nullptr;
    }
    
    // stream constructor
    Element(sbs::ifstream &ifs):
    mQuadTag(ifs.get<int>()),
    mGradQuad(std::make_unique<GradQuad>(ifs)),
    mPRT(ifs.get<bool>() ? std::make_unique<PRT>(ifs) : nullptr) {
        mTransform = nullptr;
    }
    
    // destructor
    virtual ~Element() = default;
    
    // type info
    virtual std::string typeInfo() const = 0;
    
    
    /////////////////////////// point ///////////////////////////
    // get point
    virtual Point &getPoint(int ipnt) const = 0;
    
    // get point nu
    std::array<int, spectral::nPEM> getPointNu_1() const;
    
    // point set
    void pointSet(bool elemInFourier);
    
    // get nr
    int getNr() const {
        return mNr;
    }
    
    // get nu
    int getNu() const {
        return mNu;
    }
    
    // get quad tag
    int getQuadTag() const {
        return mQuadTag;
    }
    
    // find boundary points by tag
    std::vector<int>
    findBoundaryPointsByTag(const std::vector<int> &boundaryMeshTags) const;
    
    // find boundary points by crds
    std::vector<int>
    findBoundaryPointsByCrds(const std::vector<double> &boundaryCrdsRorZ,
                             const std::vector<double> &boundaryCrdsTorS,
                             double distTol) const;
    
    /////////////////////////// domain ///////////////////////////
    // set domain tag
    void setDomainTag(int domainTag) {
        mDomainTag = domainTag;
    }
    
    // get domain tag
    int getDomainTag() const {
        return mDomainTag;
    }
    
    
protected:
    /////////////////////////// crd transform ///////////////////////////
    // set to RTZ by material or PRT
    void setToRTZ_ByMaterialOrPRT();
    
    // create coordinate transform
    // internally called by setToRTZ_ByMaterialOrPRT
    // externally called by prepareWavefieldOutput
    void createCoordTransform();
    
    
protected:
    // tag in the spectral-element mesh
    // this tag is mpi-independent and unique across processors
    const int mQuadTag;
    
    // gradient operator
    const std::unique_ptr<const GradQuad> mGradQuad;
    
    // particle relabelling
    const std::unique_ptr<const PRT> mPRT;
    
    // order
    int mNr = 0;
    int mNu = 0;
    
    // tag (position) in the computational domain
    // this tag is mpi-dependent
    int mDomainTag = -1;
    
    // coordinate transform
    std::unique_ptr<const CoordTransform> mTransform;
    
    
    ////////////////////////////////////////
    //////////////// static ////////////////
    ////////////////////////////////////////
    
public:
    // stream-based creator
    static std::shared_ptr<Element>
    createFromStream(sbs::ifstream &ifs, const Domain &domain);
};

#endif /* Element_hpp */
