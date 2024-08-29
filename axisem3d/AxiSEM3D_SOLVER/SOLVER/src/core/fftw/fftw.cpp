//
//  fftw.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/22/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  global FFTW solvers for core

#include "fftw.hpp"

// to use template member functions implemented in cpp
#include "SolverFFTW.cpp"

namespace fftw {
    // global FFTW solvers
    SolverFFTW<Real, 1> gFFT_1;
    SolverFFTW<Real, 3> gFFT_3;
    SolverFFTW<Real, nPEM * 3> gFFT_N3;
    SolverFFTW<Real, nPEM * 6> gFFT_N6;
    SolverFFTW<Real, nPEM * 9> gFFT_N9;
    
    // create plans
    void createPlans(double timeLimitForPlanning) {
        // begin
        timer::gPreloopTimer.begin("Creating FFT solvers");
        
        // split time limit
        double tf1 = gFFT_1.timeFactorForPlanning();
        double tf3 = gFFT_3.timeFactorForPlanning();
        double tfN3 = gFFT_N3.timeFactorForPlanning();
        double tfN6 = gFFT_N6.timeFactorForPlanning();
        double tfN9 = gFFT_N9.timeFactorForPlanning();
        double tfSum = std::max(tf1 + tf3 + tfN3 + tfN6 + tfN9,
                                numerical::dEpsilon);
        double timeLimitUnit = timeLimitForPlanning / tfSum;
        
        // 1
        timer::gPreloopTimer.begin("FFT with HOWMANY = 1");
        gFFT_1.createPlans(timeLimitUnit * tf1);
        gFFT_1.diagnosticInfo();
        timer::gPreloopTimer.ended("FFT with HOWMANY = 1");
        
        // 3
        timer::gPreloopTimer.begin("FFT with HOWMANY = 3");
        gFFT_3.createPlans(timeLimitUnit * tf3);
        gFFT_3.diagnosticInfo();
        timer::gPreloopTimer.ended("FFT with HOWMANY = 3");
        
        // N3
        timer::gPreloopTimer.begin("FFT with HOWMANY = # GLL per Element * 3");
        gFFT_N3.createPlans(timeLimitUnit * tfN3);
        gFFT_N3.diagnosticInfo();
        timer::gPreloopTimer.ended("FFT with HOWMANY = # GLL per Element * 3");
        
        // N6
        timer::gPreloopTimer.begin("FFT with HOWMANY = # GLL per Element * 6");
        gFFT_N6.createPlans(timeLimitUnit * tfN6);
        gFFT_N6.diagnosticInfo();
        timer::gPreloopTimer.ended("FFT with HOWMANY = # GLL per Element * 6");
        
        // N9
        timer::gPreloopTimer.begin("FFT with HOWMANY = # GLL per Element * 9");
        gFFT_N9.createPlans(timeLimitUnit * tfN9);
        gFFT_N9.diagnosticInfo();
        timer::gPreloopTimer.ended("FFT with HOWMANY = # GLL per Element * 9");
        
        // end
        timer::gPreloopTimer.ended("Creating FFT solvers");
    }
    
    // clear plans
    void clearPlans() {
        gFFT_1.clearPlans();
        gFFT_3.clearPlans();
        gFFT_N3.clearPlans();
        gFFT_N6.clearPlans();
        gFFT_N9.clearPlans();
    }
    
    // verbose
    std::string verbose() {
        using namespace bstring;
        std::stringstream ss;
        ss << boxTitle("FFT Solvers");
        ss << boxEquals(0, 43, "real number precision", typeName<Real>());
        ss << boxEquals(0, 43, "\"howmany\" on GLL-points", "{1, 3}");
        ss << boxEquals(0, 43, "\"howmany\" on elements", "{N3, N6, N9}");
        ss << gFFT_1.verbose("FFT on GLL-points, \"howmany\" = 1");
        ss << gFFT_3.verbose("FFT on GLL-points, \"howmany\" = 3");
        ss << gFFT_N3.verbose("FFT on elements, \"howmany\" = N3");
        ss << gFFT_N6.verbose("FFT on elements, \"howmany\" = N6");
        ss << gFFT_N9.verbose("FFT on elements, \"howmany\" = N9");
        ss << boxBaseline() << "\n\n";
        return ss.str();
    }
}
