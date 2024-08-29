//
//  SourceTimeFunction.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/5/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  source time function

#ifndef SourceTimeFunction_hpp
#define SourceTimeFunction_hpp

#include "eigen.hpp"
#include "spectral.hpp"
#include "vector_tools.hpp"

template <typename T, int ndim>
class SourceTimeFunction {
public:
    // pattern matrix
    typedef Eigen::Matrix<T, Eigen::Dynamic, 1> TColX;
    typedef Eigen::Matrix<std::complex<T>, Eigen::Dynamic,
    spectral::nPEM * ndim> CTMatXND;
    
    // constructor for single pattern
    SourceTimeFunction(bool alignedToTimeStep, const std::vector<double> &times,
                       const TColX &singleSTF, const CTMatXND &singlePattern):
    mAlignedToTimeStep(alignedToTimeStep), mTimePoints(times),
    mSingle(true),
    mSingleSTF(singleSTF), mSinglePattern(singlePattern),
    mMultiplePatterns() {
        checkTimePoints();
        checkPatterns();
    }
    
    // constructor for multiple patterns
    SourceTimeFunction(bool alignedToTimeStep, const std::vector<double> &times,
                       const std::vector<CTMatXND> &multiplePatterns):
    mAlignedToTimeStep(alignedToTimeStep), mTimePoints(times),
    mSingle(false),
    mSingleSTF(TColX(0, 1)), mSinglePattern(CTMatXND(0, spectral::nPEM * ndim)),
    mMultiplePatterns(multiplePatterns) {
        checkTimePoints();
        checkPatterns();
    }
    
    // get value at a time step
    void getPatternAtTimeStep(int timestep, double time,
                              CTMatXND &result) const {
        int nu_1 = getPatternNu_1();
        if (mAlignedToTimeStep) {
            // aligned
            if (mSingle) {
                result.topRows(nu_1) = mSingleSTF[timestep] * mSinglePattern;
            } else {
                result.topRows(nu_1) = mMultiplePatterns[timestep];
            }
        } else {
            // non-aligned
            int index0 = -1, index1 = -1;
            double factor0 = 0., factor1 = 0.;
            vector_tools::linearInterpSorted(mTimePoints, time, index0, index1,
                                             factor0, factor1);
            if (mSingle) {
                result.topRows(nu_1) = ((mSingleSTF[index0] * (T)factor0 +
                                         mSingleSTF[index1] * (T)factor1) *
                                        mSinglePattern);
            } else {
                result.topRows(nu_1) = (mMultiplePatterns[index0] * (T)factor0 +
                                        mMultiplePatterns[index1] * (T)factor1);
            }
        }
    }
    
    // get order
    int getPatternNu_1() const {
        if (mSingle) {
            return (int)mSinglePattern.rows();
        } else {
            return (int)mMultiplePatterns[0].rows();
        }
    }
    
private:
    // check time points
    void checkTimePoints() const {
        // use adjacent_find and greater_equal to check if
        // the time points are strictly increasing
        auto adjIter = std::adjacent_find(mTimePoints.begin(),
                                          mTimePoints.end(),
                                          std::greater_equal<double>());
        if (adjIter != mTimePoints.end()) {
            throw std::runtime_error("SourceTimeFunction::checkTimePoints || "
                                     "Time points are not increasing.");
        }
        // check time size
        if (mTimePoints.size() < 2) {
            throw std::runtime_error("SourceTimeFunction::checkTimePoints || "
                                     "Too few time points.");
        }
    }
    
    // check patterns
    void checkPatterns() const {
        if (mSingle) {
            if (mTimePoints.size() != mSingleSTF.size()) {
                throw std::runtime_error("SourceTimeFunction::checkPatterns || "
                                         "Incompatible factor counts in "
                                         "single-pattern STF.");
            }
        } else {
            if (mTimePoints.size() != mMultiplePatterns.size()) {
                throw std::runtime_error("SourceTimeFunction::checkPatterns || "
                                         "Incompatible pattern counts in "
                                         "multiple-pattern STF.");
            }
            // make sure all the patterns have the same shape
            for (const CTMatXND &pattern: mMultiplePatterns) {
                if (pattern.rows() != mMultiplePatterns[0].rows()) {
                    throw std::runtime_error("SourceTimeFunction::checkPatterns"
                                             " || Incompatible pattern shapes "
                                             "in multiple-pattern STF.");
                }
            }
        }
    }
    
private:
    ////////////// time //////////////
    // aligned to timesteps
    // true: directly retrieve stf using timestep as index
    // false: interpolation in time is needed
    const bool mAlignedToTimeStep;
    
    // time points
    const std::vector<double> mTimePoints;
    
    ////////////// pattern //////////////
    // pattern flag
    const bool mSingle;
    
    // single pattern
    const TColX mSingleSTF;
    const CTMatXND mSinglePattern;
    
    // multiple patterns
    const std::vector<CTMatXND> mMultiplePatterns;
};

#endif /* SourceTimeFunction_hpp */
