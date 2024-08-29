// Attenuation.h
// created by Kuangdai on 27-Apr-2016 
// base class of attenuation based on SLS (standard linear solids)

#pragma once

#include "eigenc.h"
#include "global.h"
#include "sbs.hpp"

class Attenuation {
public:
    
    Attenuation(int nsls, const RColX &alpha, 
        const RColX &beta, const RColX &gamma): 
        mNSLS(nsls), mAlpha(alpha), 
        mBeta(beta), mGamma(gamma) {};
    
    virtual ~Attenuation() {};
    
    // check memory variable size
    virtual void checkCompatibility(int Nr) const = 0;
    
    // reset to zero 
    virtual void resetZero() = 0; 
    
    virtual void dumpToStream(sbs::ofstream &ofs, const RRowN &ifact) const {
        ofs<<mAlpha<<mBeta<<mGamma;
    }
    
protected:
    int mNSLS;
    RColX mAlpha;
    RColX mBeta;
    RColX mGamma;        
};

