//
//  fftw.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/22/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  global FFTW solvers for core

#ifndef fftw_hpp
#define fftw_hpp

#include "SolverFFTW.hpp"
#include "spectral.hpp"
#include "numerical.hpp"

namespace fftw {
    using numerical::Real;
    using spectral::nPEM;
    
    // global FFTW solvers
    extern SolverFFTW<Real, 1> gFFT_1;
    extern SolverFFTW<Real, 3> gFFT_3;
    extern SolverFFTW<Real, nPEM * 3> gFFT_N3;
    extern SolverFFTW<Real, nPEM * 6> gFFT_N6;
    extern SolverFFTW<Real, nPEM * 9> gFFT_N9;
    
    // create plans
    void createPlans(double timeLimitForPlanning);
    
    // clear plans
    void clearPlans();
    
    // verbose
    std::string verbose();
}

#endif /* fftw_hpp */
