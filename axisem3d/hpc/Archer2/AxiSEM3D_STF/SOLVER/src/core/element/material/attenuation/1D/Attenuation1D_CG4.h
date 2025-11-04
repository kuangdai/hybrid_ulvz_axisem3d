// Attenuation1D_CG4.h
// created by Kuangdai on 29-Apr-2016 
// 1D attenuation on coarse grid

#pragma once

#include "Attenuation1D.h"
#include "eigen_cg4.h"
#include <vector>

class Attenuation1D_CG4: public Attenuation1D {    
public:
    
    Attenuation1D_CG4(int nsls, const RColX &alpha, 
        const RColX &beta, const RColX &gamma, int Nu, 
        const RRow4 &dkappa, const RRow4 &dmu, bool doKappa);
        
    // STEP 2.1: R ==> stress
    void applyToStress(vec_ar6_CMatPP &stress) const;

    // STEP 2.3: strain ==> R
    void updateMemoryVariables(const vec_ar6_CMatPP &strain);
    
    // check memory variable size
    void checkCompatibility(int Nr) const;
    
    // reset to zero 
    void resetZero(); 
    
    void dumpToStream(sbs::ofstream &ofs, const RRowN &ifact) const {
        ofs<<std::string("AttenuationCG4");
        ofs<<true;
        Attenuation::dumpToStream(ofs, ifact);
        
        int cg4[4];
        cg4[0]=1*5+1;
        cg4[1]=1*5+3;
        cg4[2]=3*5+1;
        cg4[3]=3*5+3;
        
        
        Eigen::Matrix<Real, 2,2, Eigen::RowMajor> temp;
        int ipnt;
        ipnt = 0;
        for (int ipol = 0; ipol < 2; ipol++) {
            for (int jpol = 0; jpol < 2; jpol++) {
                temp(ipol,jpol) = (mDKappa3(ipnt)/3-mDMu2(ipnt)/3) ;
                ipnt++;
            }
        }
        ofs << true<<temp<<0<<4;
        
        ipnt = 0;
        for (int ipol = 0; ipol < 2; ipol++) {
            for (int jpol = 0; jpol < 2; jpol++) {
                temp(ipol,jpol) = mDMu(ipnt) ;
                ipnt++;
            }
        }
        ofs << true<<temp<<0<<4;
        
        #ifndef _SAVE_MEMORY
        ipnt = 0;
        for (int ipol = 0; ipol < 2; ipol++) {
            for (int jpol = 0; jpol < 2; jpol++) {
                temp(ipol,jpol) = mDMu(ipnt) ;
                ipnt++;
            }
        }
        ofs << true<<temp<<0<<4;
        
        #endif
        
    }
    
    
    
private:
    
    // memory variables
    vec_ar6_CRow4 mStressR;
    std::vector<vec_ar6_CRow4> mMemVar;
    // modules
    RRow4 mDKappa3; // dkappa * 3
    RRow4 mDMu;     // dmu
    RRow4 mDMu2;    // dmu * 2
    bool mDoKappa;
};