//
//  Attenuation.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 2/25/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  attenuation

#ifndef Attenuation_hpp
#define Attenuation_hpp

// dmu * 2
#ifdef _SAVE_MEMORY
// computed on the fly as static variables
#define xDMu2 sDMu2
#else
// precomputed and stored as member variables
#define xDMu2 mDMu2
#endif

#include "eigen_element.hpp"
#include "FieldArithmetic.hpp"

namespace eigen {
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> DColX;
    typedef Eigen::Matrix<numerical::Real, Eigen::Dynamic, 1> RColX;
}

class Attenuation {
public:
    // parameter constructor
    Attenuation(bool is1D,
                const eigen::DColX &alpha,
                const eigen::DColX &betta,
                const eigen::DColX &gamma):
    m1D(is1D),
    mAlpha(alpha.cast<numerical::Real>()),
    mBetta(betta.cast<numerical::Real>()),
    mGamma(gamma.cast<numerical::Real>()) {
        // nothing
    }
    
    // copy constructor
    Attenuation(const Attenuation &other):
    m1D(other.m1D),
    mAlpha(other.mAlpha),
    mBetta(other.mBetta),
    mGamma(other.mGamma) {
        // nothing
    }
    
    // stream constructor
    Attenuation(sbs::ifstream &ifs):
    m1D(ifs.get<bool>()),
    mAlpha(ifs.get<eigen::RColX>()),
    mBetta(ifs.get<eigen::RColX>()),
    mGamma(ifs.get<eigen::RColX>()) {
        // nothing
    }
    
    // destructor
    virtual ~Attenuation() = default;
    
    // clone for copy constructor
    virtual std::unique_ptr<Attenuation> clone() const = 0;
    
    // get number of SLS
    int nsls() const {
        return (int)mAlpha.rows();
    }
    
    // check compatibility
    virtual void
    checkCompatibility(int nr, bool elemInFourier, bool elastic1D) {
        // 1D/3D
        if (elastic1D != m1D) {
            throw std::runtime_error("AttenuationFull::checkCompatibility || "
                                     "Incompatible 1D/3D flags.");
        }
    }
    
    
    //////////////////////// apply ////////////////////////
    // apply attenuation in Fourier space
    virtual void apply(const eigen::vec_ar6_CMatPP_RM &strain,
                       eigen::vec_ar6_CMatPP_RM &stress, int nu_1) = 0;
    
    // apply attenuation in cardinal space
    virtual void apply(const eigen::RMatXN6 &strain,
                       eigen::RMatXN6 &stress, int nr) = 0;
    
    
    //////////////////////// data ////////////////////////
protected:
    // 1D/3D flag
    const bool m1D;
    
private:
    // coeffs
    const eigen::RColX mAlpha;
    const eigen::RColX mBetta;
    const eigen::RColX mGamma;
    
    
    /////////////////////////////////////////////////////////////////////
    /////////////////// template-based implementation ///////////////////
    /////////////////////////////////////////////////////////////////////
    
protected:
    ///////////////// update memory variable, phase 1 /////////////////
    // CASE == _1D_FR
    template <CaseFA CASE, class FMat6>
    typename std::enable_if<CASE == CaseFA::_1D_FR, void>::type
    updateMemAlphaBeta(const FMat6 &dStress, std::vector<FMat6> &memVar,
                       int nx) const {
        for (int isls = 0; isls < nsls(); isls++) {
            for (int idim = 0; idim < 6; idim++) {
                memVar[isls][nx][idim] =
                mAlpha[isls] * memVar[isls][nx][idim] +
                mBetta[isls] * dStress[nx][idim];
            }
        }
    }
    // CASE != _1D_FR
    template <CaseFA CASE, class FMat6>
    typename std::enable_if<CASE != CaseFA::_1D_FR, void>::type
    updateMemAlphaBeta(const FMat6 &dStress, std::vector<FMat6> &memVar,
                       int nx) const {
        for (int isls = 0; isls < nsls(); isls++) {
            memVar[isls] = (mAlpha[isls] * memVar[isls] +
                            mBetta[isls] * dStress);
        }
    }
    
    
    ///////////////// update memory variable, phase 2 /////////////////
    // CASE == _1D_FR
    template <CaseFA CASE, class FMat6>
    typename std::enable_if<CASE == CaseFA::_1D_FR, void>::type
    updateMemGamma(const FMat6 &dStress, std::vector<FMat6> &memVar,
                   int nx) const {
        for (int isls = 0; isls < nsls(); isls++) {
            for (int idim = 0; idim < 6; idim++) {
                memVar[isls][nx][idim] += mGamma[isls] * dStress[nx][idim];
            }
        }
    }
    // CASE != _1D_FR
    template <CaseFA CASE, class FMat6>
    typename std::enable_if<CASE != CaseFA::_1D_FR, void>::type
    updateMemGamma(const FMat6 &dStress, std::vector<FMat6> &memVar,
                   int nx) const {
        for (int isls = 0; isls < nsls(); isls++) {
            memVar[isls] += mGamma[isls] * dStress;
        }
    }
    
    
    ////////////////////////////////////////
    //////////////// static ////////////////
    ////////////////////////////////////////
    
public:
    // stream-based creator
    static std::unique_ptr<Attenuation> createFromStream(sbs::ifstream &ifs);
};

#endif /* Attenuation_hpp */
