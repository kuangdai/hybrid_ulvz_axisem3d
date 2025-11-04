// SFCoupling3D.h
// created by Kuangdai on 5-Apr-2016 
// 3D solid-fluid boundary condition

#pragma once

#include "SFCoupling.h"

class SFCoupling3D: public SFCoupling {
public:   
    SFCoupling3D(const RMatX3 &n, const RMatX3 &n_invmf): 
        mNormal_unassembled(n),
        mNormal_assembled_invMassFluid(n_invmf) {};
    
    // solid-fluid coupling
    void coupleFluidToSolid(const CColX &fluidStiff, CMatX3 &solidStiff) const; 
    void coupleSolidToFluid(const CMatX3 &solidDispl, CColX &fluidStiff) const;
    
    // verbose
    std::string verbose() const {return "SFCoupling3D";};
    
    // check compatibility
    void checkCompatibility(int nr) const;
    
    void dumpToStream(sbs::ofstream &ofs, SolidPoint *sp, FluidPoint *fp) const{
        ofs<<"SolidFluidCoupling3D";
        ofs<<sp->mSolidTag;
        ofs<<fp->mFluidTag;
        ofs<<mNormal_unassembled<<mNormal_assembled_invMassFluid;
    }  

private:    
    RMatX3 mNormal_unassembled;
    RMatX3 mNormal_assembled_invMassFluid;
};
