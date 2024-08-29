// SolverFFTW_1.h
// created by Kuangdai on 21-Apr-2016 
// perform FFT using fftw

#pragma once

#include "SolverFFTW.h"

class SolverFFTW_1_nonstatic {
public:
    // initialize plans
    void initialize(int Nmax);
    // finalize plans
    void finalize();
    
    // get input and output
    RColX &getR2C_RMat() {return sR2C_RMat;};
    CColX &getR2C_CMat() {return sR2C_CMat;};
    RColX &getC2R_RMat() {return sC2R_RMat;};    
    CColX &getC2R_CMat() {return sC2R_CMat;};
     
    // forward, real => complex
    void computeR2C(int nr);
    // backward, complex => real
    void computeC2R(int nr);
        
private:
    int sNmax = 0;
    std::vector<PlanFFTW> sR2CPlans;
    std::vector<PlanFFTW> sC2RPlans;
    RColX sR2C_RMat;
    CColX sR2C_CMat;
    RColX sC2R_RMat;
    CColX sC2R_CMat;
};
