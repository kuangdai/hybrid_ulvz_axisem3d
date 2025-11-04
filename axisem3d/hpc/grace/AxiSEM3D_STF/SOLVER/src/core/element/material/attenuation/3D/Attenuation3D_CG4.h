// Attenuation3D_CG4.h
// created by Kuangdai on 29-Apr-2016 
// 3D attenuation on coarse grid 

#pragma once

#include "Attenuation3D.h"
#include "eigen_cg4.h"
#include <vector>

class Attenuation3D_CG4: public Attenuation3D {
public:
    
    Attenuation3D_CG4(int nsls, 
        const RColX &alpha, const RColX &beta, const RColX &gamma, 
        const RMatX4 &dkappa, const RMatX4 &dmu, bool doKappa);
    
    // STEP 2.1: R ==> stress
    void applyToStress(RMatXN6 &stress) const;

    // STEP 2.3: strain ==> R
    void updateMemoryVariables(const RMatXN6 &strain);
    
    // check memory variable size
    void checkCompatibility(int Nr) const;
    
    // reset to zero 
    void resetZero(); 
    
    void dumpToStream(sbs::ofstream &ofs, const RRowN &ifact) const {
        ofs<<std::string("AttenuationCG4");
        ofs<<false;
        Attenuation::dumpToStream(ofs, ifact);
        
        Eigen::Matrix<Real, 1, 4> if4;
        if4(0) = ifact(1*5+1);
        if4(1) = ifact(1*5+3);
        if4(2) = ifact(3*5+1);
        if4(3) = ifact(3*5+3);
        
        RMatX4 temp = ((mDKappa3 - mDMu2) /3);
        // temp.array().rowwise() /= if4.array();
        ofs << false<<temp;
        
        temp = mDMu;
        // temp.array().rowwise() /= if4.array();
        ofs << false<<temp;
        
        #ifndef _SAVE_MEMORY
        temp = mDMu2;
        // temp.array().rowwise() /= if4.array();
        ofs << false<<temp;
        
        #endif
        
    }
    
private:
    
    // memory variables
    RMatX46 mStressR;
    RMatX46 mStrain4;
    std::vector<RMatX46> mMemVar;
    // modules
    RMatX4 mDKappa3; // dkappa * 3
    RMatX4 mDMu;     // dmu
    RMatX4 mDMu2;     // dmu
    bool mDoKappa;
};
