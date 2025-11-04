//
//  Point.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/24/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  GLL point

#ifndef Point_hpp
#define Point_hpp

#include "eigen_point.hpp"
#include "Mass.hpp"

class TimeScheme;

class Point {
public:
    // parameter constructor
    Point(int nr, const eigen::DRow2 &crds, int meshTag,
          std::unique_ptr<const Mass> &mass):
    mNr(nr), mNu(mNr / 2), mCoords(crds),
    mMeshTag(meshTag), mMass(mass.release()) {
        // nothing
    }
    
    // stream constructor
    Point(sbs::ifstream &ifs):
    mNr(ifs.get<int>()), mNu(mNr / 2), mCoords(ifs.get<eigen::DRow2>()),
    mMeshTag(ifs.get<int>()), mMass(Mass::createFromStream(ifs)) {
        // nothing
    }
    
    // destructor
    virtual ~Point() = default;
    
    // type info
    std::string typeInfo() const;
    
    
    /////////////////////////// properties ///////////////////////////
    // nr
    int getNr() const {
        return mNr;
    }
    
    // nu
    int getNu() const {
        return mNu;
    }
    
    // location
    const eigen::DRow2 &getCoords() const {
        return mCoords;
    }
    
    // tag
    int getMeshTag() const {
        return mMeshTag;
    }
    
    
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
    // order
    const int mNr;
    const int mNu;
    
    // location (s, z)
    const eigen::DRow2 mCoords;
    
    // tag in the spectral-element mesh, determined by location
    // * a pair of solid and fluid points at the same location
    //   have the same mesh tag
    // TODO: make this tag mpi-independent
    const int mMeshTag;
    
    // mass
    const std::unique_ptr<const Mass> mMass;
    
    // tag (index) in the computational domain, mpi-dependent
    int mDomainTag = -1;
    
    
    ////////////////////////////////////////
    //////////////// static ////////////////
    ////////////////////////////////////////
    
public:
    // stream-based creator
    static std::shared_ptr<Point>
    createFromStream(sbs::ifstream &ifs, const TimeScheme &timeScheme);
};

#endif /* Point_hpp */
