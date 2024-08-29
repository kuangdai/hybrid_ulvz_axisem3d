//
//  numerical.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/10/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  numerical precision and constants

#ifndef numerical_hpp
#define numerical_hpp

#include <complex>
#include <limits>

namespace numerical {
    // large int
    typedef std::ptrdiff_t Int;
    
    // solver precision
#ifdef _USE_DOUBLE
    typedef double Real;
#else
    typedef float Real;
#endif
    
    // epsilon
    const Real rEpsilon = std::numeric_limits<Real>::epsilon() * (Real)2.;
    const double dEpsilon = std::numeric_limits<double>::epsilon() * 2.;
    
    // complex
    typedef std::complex<Real> ComplexR;
    typedef std::complex<double> ComplexD;
    
    // constants
    const double dPi = 3.141592653589793238463;
    const double dDegree = dPi / 180.;
    
    // error or undefined number
    const double dErr = -12345e99;
}

#endif /* numerical_hpp */
