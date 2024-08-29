// Attenuation1D_Full.h
// created by Kuangdai on 29-Apr-2016 
// 1D attenuation on full grid

#pragma once

#include "Attenuation1D.h"
#include <vector>

class Attenuation1D_Full: public Attenuation1D {
public:
    
    Attenuation1D_Full(int nsls, const RColX &alpha, 
        const RColX &beta, const RColX &gamma, int Nu, 
        const RMatPP &dkappa, const RMatPP &dmu, bool doKappa);
        
    // STEP 2.1: R ==> stress
    void applyToStress(vec_ar6_CMatPP &stress) const;

    // STEP 2.3: strain ==> R
    void updateMemoryVariables(const vec_ar6_CMatPP &strain);
    
    // check memory variable size
    void checkCompatibility(int Nr) const;
    
    // reset to zero 
    void resetZero(); 
    
    void dumpToStream(sbs::ofstream &ofs, const RRowN &ifact) const {
        ofs<<std::string("AttenuationFull");
        ofs<<true;
        Attenuation::dumpToStream(ofs, ifact);
        
        RMatPP temp;
        int ipnt;
        ipnt = 0;
        for (int ipol = 0; ipol <= nPol; ipol++) {
            for (int jpol = 0; jpol <= nPol; jpol++) {
                temp(ipol,jpol) = (mDKappa3(ipol, jpol)/3-mDMu2(ipol, jpol)/3) ;
                ipnt++;
            }
        }
        ofs << true<<temp<<0<<25;
        
        ipnt = 0;
        for (int ipol = 0; ipol <= nPol; ipol++) {
            for (int jpol = 0; jpol <= nPol; jpol++) {
                temp(ipol,jpol) = mDMu(ipol, jpol) ;
                ipnt++;
            }
        }
        ofs << true<<temp<<0<<25;
        
        #ifndef _SAVE_MEMORY
        ipnt = 0;
        for (int ipol = 0; ipol <= nPol; ipol++) {
            for (int jpol = 0; jpol <= nPol; jpol++) {
                temp(ipol,jpol) = mDMu2(ipol, jpol) ;
                ipnt++;
            }
        }
        ofs << true<<temp<<0<<25;
        
        #endif
        
    }
    
private:
    // memory variables
    vec_ar6_CMatPP mStressR;
    std::vector<vec_ar6_CMatPP> mMemVar;
    // modules
    RMatPP mDKappa3; // dkappa * 3
    RMatPP mDMu;     // dmu
    RMatPP mDMu2;    // dmu * 2
    bool mDoKappa;
};

