//
//  ClaytonFluid1D.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/29/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  Clayton-Enquist ABC for fluid points in 1D

#ifndef ClaytonFluid1D_hpp
#define ClaytonFluid1D_hpp

#include "ClaytonFluid.hpp"
#include "numerical.hpp"

class ClaytonFluid1D: public ClaytonFluid {
public:
    // parameter constructor
    ClaytonFluid1D(const std::shared_ptr<FluidPoint> &fp,
                   double rho, double vp, double area);
    
    // stream constructor
    ClaytonFluid1D(sbs::ifstream &ifs, const Domain &domain);
    
    // apply ABC
    void apply() const;
    
private:
    // area / (rho * vp)
    const numerical::Real mAreaOverRhoVp;
};

#endif /* ClaytonFluid1D_hpp */
