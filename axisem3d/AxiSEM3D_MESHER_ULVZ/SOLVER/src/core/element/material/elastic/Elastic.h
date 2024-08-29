// Elastic.h
// created by Kuangdai on 2-Apr-2016 
// base class of elasticity

#pragma once

class SolidResponse;
#include <string>
#include "sbs.hpp"
#include "eigenc.h"
class Elastic {
public:
    
    virtual ~Elastic() {};
    
    // STEP 2: strain ==> stress
    virtual void strainToStress(SolidResponse &response) const = 0;
        
    // check compatibility
    virtual void checkCompatibility(int Nr) const = 0; 
    
    // verbose
    virtual std::string verbose() const = 0;
    
    // 1D or Fourier space
    virtual bool is1D() const = 0;
    
    // need TIso
    virtual bool needTIso() const = 0;
    
    // reset to zero 
    virtual void resetZero() = 0; 
    
    virtual void dumpToStream(sbs::ofstream &ofs, const RRowN &ifact) const = 0;
};
