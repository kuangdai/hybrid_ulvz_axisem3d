//
//  SolverFFTW.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/22/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  FFT solver based on FFTW3
//  v2.0 features:
//  * all logical sizes share the same static memory
//  * external data feeding and retrieving are enforced for readability
//  * template-based, available for both core and preloop

#include "SolverFFTW.hpp"
#include "vector_tools.hpp"
#include "bstring.hpp"
#include "mpi.hpp"
#include "timer.hpp"

// create plans after adding ALL NRs
template <typename floatT, int HOWMANY>
void SolverFFTW<floatT, HOWMANY>::createPlans(double timeLimitForPlanning) {
    // clear if exists
    clearPlans();
    
    // empty
    if (mToBeNRs.size() == 0) {
        return;
    }
    
    // sort and unique NRs
    vector_tools::sortUnique(mToBeNRs);
    
    // max size
    int maxNR = *(mToBeNRs.end() - 1);
    int maxNC = maxNR / 2 + 1;
    
    // data allocation
    mRMatXM.resize(maxNR, HOWMANY);
    mCMatXM.resize(maxNC, HOWMANY);
    
    // split time limit by NR
    double factor = 0.;
    for (const int &NR: mToBeNRs) {
        // use NR * log2(NR) for NR scaling
        factor += NR * std::log2(NR * 1.);
    }
    // factor * 2: R2C and C2R
    double timeUnit = timeLimitForPlanning / (factor * 2.);
    
    // index mapping
    mPlanIndex = std::vector<int>(maxNR + 1, -1);
    
    // create plans
    int loc = 0;
    mPlansR2C.reserve(mToBeNRs.size());
    mPlansC2R.reserve(mToBeNRs.size());
    for (const int &NR: mToBeNRs) {
        // split time limit by NR
        InterfaceFFTW<floatT>::timelimit(timeUnit * NR * std::log2(NR * 1.));
        
        // r2c plan
        typename InterfaceFFTW<floatT>::Plan plan_r2c =
        InterfaceFFTW<floatT>::planR2C(1, &NR, HOWMANY,
                                       mRMatXM.data(), NULL, 1, maxNR,
                                       reinterpret_cast<typename
                                       InterfaceFFTW<floatT>::fftw_cmplx *>
                                       (mCMatXM.data()), NULL, 1, maxNC,
                                       FFTW_PATIENT);
        mPlansR2C.push_back(plan_r2c);
        
        // c2r plan
        typename InterfaceFFTW<floatT>::Plan plan_c2r =
        InterfaceFFTW<floatT>::planC2R(1, &NR, HOWMANY,
                                       reinterpret_cast<typename
                                       InterfaceFFTW<floatT>::fftw_cmplx *>
                                       (mCMatXM.data()), NULL, 1, maxNC,
                                       mRMatXM.data(), NULL, 1, maxNR,
                                       FFTW_PATIENT);
        mPlansC2R.push_back(plan_c2r);
        
        // index mapping
        mPlanIndex[NR] = loc++;
    }
    
    // NR list
    mCurrentNRs = mToBeNRs;
    mToBeNRs.clear();
}

// clear plans
template <typename floatT, int HOWMANY>
void SolverFFTW<floatT, HOWMANY>::clearPlans() {
    // empty
    if (mCurrentNRs.size() == 0) {
        return;
    }
    
    // index mapping
    mPlanIndex = std::vector<int>(1, -1);
    
    // destroy plans
    for (int iplan = 0; iplan < mPlansR2C.size(); iplan++) {
        InterfaceFFTW<floatT>::destroy(mPlansR2C[iplan]);
        InterfaceFFTW<floatT>::destroy(mPlansC2R[iplan]);
    }
    mPlansR2C.clear();
    mPlansC2R.clear();
    
    // data
    mRMatXM.resize(0, HOWMANY);
    mCMatXM.resize(0, HOWMANY);
    
    // NR list
    mCurrentNRs.clear();
}

// time factor for planning
template <typename floatT, int HOWMANY>
double SolverFFTW<floatT, HOWMANY>::timeFactorForPlanning() const {
    // empty
    if (mToBeNRs.size() == 0) {
        return 0.;
    }
    
    // mToBeNRs is unsorted here
    std::vector<int> mtempNRs(mToBeNRs);
    vector_tools::sortUnique(mtempNRs);
    
    // accumulate time factor
    double factor = 0.;
    for (const int &NR: mtempNRs) {
        // use NR * log2(NR) for NR scaling
        factor += NR * std::log2(NR * 1.);
    }
    
    // use HOWMANY ^ 0.75 for HOWMANY scaling
    return factor * std::pow(HOWMANY * 1., 0.75);
}

// diagnostic info after creating plans
template <typename floatT, int HOWMANY>
void SolverFFTW<floatT, HOWMANY>::diagnosticInfo() const {
    // statistics
    int numNR_MaxMPI = 0;
    int maxNR_MaxMPI = 0;
    int sumNR_MaxMPI = 0;
    double memGB_MaxMPI = 0.;
    statistics(numNR_MaxMPI, maxNR_MaxMPI, sumNR_MaxMPI, memGB_MaxMPI);
    
    // write diagnostic info
    timer::gPreloopTimer.message("HOWMANY = " + bstring::toString(HOWMANY));
    timer::gPreloopTimer.message("# NR's (max over MPI ranks) = " +
                                 bstring::toString(numNR_MaxMPI));
    timer::gPreloopTimer.message("max NR (max over MPI ranks) = " +
                                 bstring::toString(maxNR_MaxMPI));
    timer::gPreloopTimer.message("sum NR (max over MPI ranks) = " +
                                 bstring::toString(sumNR_MaxMPI));
    timer::gPreloopTimer.message("*** non-scalable memory allocation "
                                 "(max over MPI ranks) = " +
                                 bstring::toString(memGB_MaxMPI) + " GB ***");
}

// verbose
template <typename floatT, int HOWMANY> std::string
SolverFFTW<floatT, HOWMANY>::verbose(const std::string &subtitle) const {
    // statistics
    int numNR_MaxMPI = 0;
    int maxNR_MaxMPI = 0;
    int sumNR_MaxMPI = 0;
    double memGB_MaxMPI = 0.;
    statistics(numNR_MaxMPI, maxNR_MaxMPI, sumNR_MaxMPI, memGB_MaxMPI);
    
    // verbose
    std::stringstream ss;
    ss << bstring::boxSubTitle(0, subtitle);
    ss << bstring::boxEquals(2, 41, "\"howmany\"", HOWMANY);
    ss << bstring::boxEquals(2, 41, "# logical sizes (NR) "
                             "(max over MPI ranks)", numNR_MaxMPI);
    ss << bstring::boxEquals(2, 41, "maximum logical size "
                             "(max over MPI ranks)", maxNR_MaxMPI);
    ss << bstring::boxEquals(2, 41, "sum of logical sizes "
                             "(max over MPI ranks)", sumNR_MaxMPI);
    return ss.str();
}

// statistics for diagnosis and verbose
template <typename floatT, int HOWMANY>
void SolverFFTW<floatT, HOWMANY>::statistics(int &numNR_MaxMPI,
                                             int &maxNR_MaxMPI,
                                             int &sumNR_MaxMPI,
                                             double &memGB_MaxMPI) const {
    // local
    int numNR = (int)mCurrentNRs.size();
    int maxNR = (numNR == 0) ? 0 : *(mCurrentNRs.end() - 1);
    int sumNR = vector_tools::sum(mCurrentNRs);
    double memGB = sizeof(floatT) * maxNR * HOWMANY * 2. / 1e9;
    
    // mpi global
    numNR_MaxMPI = mpi::max(numNR);
    maxNR_MaxMPI = mpi::max(maxNR);
    sumNR_MaxMPI = mpi::max(sumNR);
    memGB_MaxMPI = mpi::max(memGB);
}
