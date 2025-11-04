// MassOcean1D.h
// created by Kuangdai on 3-Apr-2016 
// 1D mass with ocean

#pragma once

#include "Mass.h"

class MassOcean1D : public Mass {
public:
    MassOcean1D(double mass, double massOcean, double theta);
    
    // compute accel in-place
    void computeAccel(CMatX3 &stiff) const;
    void computeAccel(CColX &stiff) const;
    
    // verbose
    std::string verbose() const {return "MassOcean1D";};
    
    // stream dumper
    void dumpToStream(sbs::ofstream &ofs) const {
        ofs << std::string("MassOceanLoad1D");
        // IMH Cos[t]^2 + IMV Sin[t]^2
        const Real mIMH_CosT2_p_IMV_SinT2 = mInvMassR * mCost*mCost+mInvMassZ * mSint*mSint;
        // IMH Sin[t]^2 + IMV Cos[t]^2
        const Real mIMH_SinT2_p_IMV_CosT2 = mInvMassR * mSint*mSint+mInvMassZ * mCost*mCost;
        // (IMV - IMH) Cos[t] Sin[t]
        const Real mIMV_m_IMH_x_CosT_SinT = (mInvMassZ - mInvMassR) * mCost * mSint;
        // IMH (for transverse component)
        const Real mIMH = mInvMassR;
        ofs << mIMH_CosT2_p_IMV_SinT2<<mIMH_SinT2_p_IMV_CosT2<<mIMV_m_IMH_x_CosT_SinT<<mIMH;
    }
    
private:
    // the scalar mass
    Real mInvMassZ;
    Real mInvMassR;
    Real mSint;
    Real mCost;
};
