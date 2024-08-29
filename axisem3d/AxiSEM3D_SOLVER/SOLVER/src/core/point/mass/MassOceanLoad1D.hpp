//
//  MassOceanLoad1D.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/25/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  1D mass with ocean load
//  theta: angle between surface normal and z-axis

#ifndef MassOceanLoad1D_hpp
#define MassOceanLoad1D_hpp

#include "Mass.hpp"

class MassOceanLoad1D: public Mass {
public:
    // parameter constructor
    MassOceanLoad1D(double mass, double massOcean, double theta);
    
    // stream constructor
    MassOceanLoad1D(sbs::ifstream &ifs);
    
    // check compatibility
    void checkCompatibility(int nr, bool solid) const;
    
    // compute accel in-place for fluid
    void computeAccel(eigen::CColX &stiff1) const {
        throw std::runtime_error("MassOceanLoad1D::computeAccel || "
                                 "Incompatible types: "
                                 "ocean load on fluid point.");
    }
    
    // compute accel in-place for solid
    void computeAccel(eigen::CMatX3 &stiff3) const;
    
private:
    // IMH = 1 / (mass)
    // IMV = 1 / (mass + massOcean)
    // IMH Cos[t]^2 + IMV Sin[t]^2
    const numerical::Real mIMH_CosT2_p_IMV_SinT2;
    // IMH Sin[t]^2 + IMV Cos[t]^2
    const numerical::Real mIMH_SinT2_p_IMV_CosT2;
    // (IMV - IMH) Cos[t] Sin[t]
    const numerical::Real mIMV_m_IMH_x_CosT_SinT;
    // IMH (for transverse component)
    const numerical::Real mIMH;
    
    
    ////////////////////////////////////////
    //////////////// static ////////////////
    ////////////////////////////////////////
    
private:
    // workspace
    inline static eigen::CColX sStiff3_col0 = eigen::CColX(0);
};

#endif /* MassOceanLoad1D_hpp */
