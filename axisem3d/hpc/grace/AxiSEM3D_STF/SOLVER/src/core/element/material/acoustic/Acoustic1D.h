// Acoustic1D.h
// created by Kuangdai on 23-Apr-2016 
// 1D acoustic

#pragma once

#include "Acoustic.h"
#include "eigenc.h"

class Acoustic1D: public Acoustic {
public:
    // constructor
    Acoustic1D(const RMatPP &KFluid): mKStruct(KFluid) {};
    
    // STEP 2: strain ==> stress
    void strainToStress(FluidResponse &response) const;
    
    // verbose
    std::string verbose() const {return "Acoustic1D";};
    
    // 1D or Fourier space
    bool is1D() const {return true;};
    
    void dumpToStream(sbs::ofstream &ofs, const RRowN &ifact) const {
        ofs<<true;
        RMatPP temp;
        int ipnt;
        ipnt = 0;
        for (int ipol = 0; ipol <= nPol; ipol++) {
            for (int jpol = 0; jpol <= nPol; jpol++) {
                temp(ipol,jpol) = mKStruct(ipol, jpol) ;
                ipnt++;
            }
        }
        ofs << true<<temp<<0<<25;
        
    }
    
private:
    RMatPP mKStruct; 
};
