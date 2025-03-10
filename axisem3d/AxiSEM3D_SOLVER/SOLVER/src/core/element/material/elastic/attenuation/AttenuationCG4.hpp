//
//  AttenuationCG4.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/13/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  attenuation on CG4 GLL points

#ifndef AttenuationCG4_hpp
#define AttenuationCG4_hpp

#include "Attenuation.hpp"

namespace eigen {
    // Fourier
    typedef Eigen::Matrix<double, 2, 2, Eigen::RowMajor> DMat22_RM;
    typedef Eigen::Matrix<numerical::Real, 2, 2, Eigen::RowMajor> RMat22_RM;
    typedef Eigen::Matrix<numerical::ComplexR, 2, 2, Eigen::RowMajor> CMat22_RM;
    typedef std::vector<std::array<CMat22_RM, 6>> vec_ar6_CMat22_RM;
    
    // cardinal
    typedef Eigen::Matrix<double, Eigen::Dynamic, 4> DMatX4;
    typedef Eigen::Matrix<numerical::Real, Eigen::Dynamic, 4> RMatX4;
    typedef Eigen::Matrix<numerical::Real, Eigen::Dynamic, 4 * 6> RMatX46;
}

namespace fa4 {
    // property class
    typedef Property<2> Property4;
    
    // operations on all GLL points
    typedef FieldArithmetic<2> FieldArithmetic4;
}

class AttenuationCG4: public Attenuation {
public:
    // 1D constructor
    AttenuationCG4(const eigen::DColX &alpha,
                   const eigen::DColX &betta,
                   const eigen::DColX &gamma,
                   const eigen::DMat22_RM &dLambda,
                   const eigen::DMat22_RM &dMu):
    Attenuation(true, alpha, betta, gamma),
    mDLambda(dLambda), mDMu(dMu)
#ifndef _SAVE_MEMORY
    , mDMu2(eigen::DMat22_RM(dMu * 2.))
#endif
    {
        // nothing
    }
    
    // 3D constructor
    AttenuationCG4(const eigen::DColX &alpha,
                   const eigen::DColX &betta,
                   const eigen::DColX &gamma,
                   const eigen::DMatX4 &dLambda,
                   const eigen::DMatX4 &dMu):
    Attenuation(false, alpha, betta, gamma),
    mDLambda(dLambda), mDMu(dMu)
#ifndef _SAVE_MEMORY
    , mDMu2(eigen::DMatX4(dMu * 2.))
#endif
    {
        // nothing
    }
    
    // copy constructor
    AttenuationCG4(const AttenuationCG4 &other):
    Attenuation(other),
    mDLambda(other.mDLambda), mDMu(other.mDMu)
#ifndef _SAVE_MEMORY
    , mDMu2(other.mDMu2)
#endif
    {
        // nothing
    }
    
    // stream constructor
    AttenuationCG4(sbs::ifstream &ifs):
    Attenuation(ifs),
    mDLambda(ifs), mDMu(ifs)
#ifndef _SAVE_MEMORY
    , mDMu2(ifs)
#endif
    {
        // nothing
    }
    
    // clone for copy constructor
    std::unique_ptr<Attenuation> clone() const {
        return std::make_unique<AttenuationCG4>(*this);
    }
    
    // check compatibility
    void checkCompatibility(int nr, bool elemInFourier, bool elastic1D) {
        // base
        Attenuation::checkCompatibility(nr, elemInFourier, elastic1D);
        
        // nPol
        if (spectral::nPol != 4) {
            throw std::runtime_error("AttenuationCG4::checkCompatibility || "
                                     "Polynomial order (nPol) is not 4.");
        }
        
        // memory variables
        int nu_1 = nr / 2 + 1;
        if (elemInFourier) {
            mDStress_FR.resize(nu_1);
            for (int alpha = 0; alpha < nu_1; alpha++) {
                mDStress_FR[alpha].fill(eigen::CMat22_RM::Zero());
            }
            mDStress_CD.resize(0, 4 * 6);
        } else {
            mDStress_FR.clear();
            mDStress_FR.shrink_to_fit();
            mDStress_CD = eigen::RMatX46::Zero(nr, 4 * 6);
        }
        mMemVar_FR = std::vector<eigen::vec_ar6_CMat22_RM>(nsls(), mDStress_FR);
        mMemVar_CD = std::vector<eigen::RMatX46>(nsls(), mDStress_CD);
        
        // parameters
        mDMu.checkCompatibility(m1D, nr, elemInFourier,
                                "mDMu", "AttenuationCG4");
        
        // workspace
        if (elemInFourier) {
            if (sStrain4_FR.size() < nu_1) {
                sStrain4_FR.resize(nu_1);
            }
        } else {
            if (sStrain4_CD.rows() < nr) {
                sStrain4_CD.resize(nr, 4 * 6);
            }
        }
#ifdef _SAVE_MEMORY
        sDMu2.expandSize(m1D, nr);
#endif
    }
    
    
    //////////////////////// apply ////////////////////////
    // apply attenuation in Fourier space
    void apply(const eigen::vec_ar6_CMatPP_RM &strain,
               eigen::vec_ar6_CMatPP_RM &stress, int nu_1) {
#ifdef _SAVE_MEMORY
        computeDMu2();
#endif
        for (int alpha = 0; alpha < nu_1; alpha++) {
            apply<CaseFA::_1D_FR>(strain, stress, alpha,
                                  mDStress_FR, mMemVar_FR, sStrain4_FR,
                                  mDLambda, mDMu, xDMu2);
        }
    }
    
    // apply attenuation in cardinal space
    void apply(const eigen::RMatXN6 &strain,
               eigen::RMatXN6 &stress, int nr) {
#ifdef _SAVE_MEMORY
        computeDMu2();
#endif
        if (m1D) {
            apply<CaseFA::_1D_CD>(strain, stress, nr,
                                  mDStress_CD, mMemVar_CD, sStrain4_CD,
                                  mDLambda, mDMu, xDMu2);
        } else {
            apply<CaseFA::_3D_CD>(strain, stress, nr,
                                  mDStress_CD, mMemVar_CD, sStrain4_CD,
                                  mDLambda, mDMu, xDMu2);
        }
    }
    
    
    ///////////////////////// properties /////////////////////////
private:
    // Cijkl
    const fa4::Property4 mDLambda;
    const fa4::Property4 mDMu;
    
    // dmu * 2
#ifdef _SAVE_MEMORY
    // computed on the fly as static variables
    inline static fa4::Property4 sDMu2;
    void computeDMu2() const {
        static const numerical::Real two = 2.;
        sDMu2.set(m1D, mDMu, two);
    }
#else
    // precomputed and stored as member variables
    const fa4::Property4 mDMu2;
#endif
    
    // memory variables in Fourier space
    eigen::vec_ar6_CMat22_RM mDStress_FR;
    std::vector<eigen::vec_ar6_CMat22_RM> mMemVar_FR;
    
    // memory variables in cardinal space
    eigen::RMatX46 mDStress_CD;
    std::vector<eigen::RMatX46> mMemVar_CD;
    
    // workspace
    eigen::vec_ar6_CMat22_RM sStrain4_FR;
    eigen::RMatX46 sStrain4_CD;
    
    
    /////////////////////////////////////////////////////////////////////
    /////////////////// template-based implementation ///////////////////
    /////////////////////////////////////////////////////////////////////
    
    /////////////////// stress -= mem ///////////////////
    // CASE == _1D_FR
    void subtractMemFromStress(eigen::vec_ar6_CMatPP_RM &stressN,
                               int nx) const {
        for (int isls = 0; isls < nsls(); isls++) {
            const eigen::vec_ar6_CMat22_RM &memVar4 = mMemVar_FR[isls];
            for (int idim = 0; idim < 6; idim++) {
                stressN[nx][idim](1, 1) -= memVar4[nx][idim](0, 0);
                stressN[nx][idim](1, 3) -= memVar4[nx][idim](0, 1);
                stressN[nx][idim](3, 1) -= memVar4[nx][idim](1, 0);
                stressN[nx][idim](3, 3) -= memVar4[nx][idim](1, 1);
            }
        }
    }
    // CASE != _1D_FR
    void subtractMemFromStress(eigen::RMatXN6 &stressN, int nx) const {
        using spectral::nPEM;
        using spectral::nPED;
        
        // seqN
        using Eigen::seqN;
        using Eigen::fix;
        const auto &seqN_Row = seqN(fix<0>, nx);
        const auto &seqN_ColN0 = seqN(fix<nPED * 1 + 1>, fix<6>, fix<nPEM>);
        const auto &seqN_ColN1 = seqN(fix<nPED * 1 + 3>, fix<6>, fix<nPEM>);
        const auto &seqN_ColN2 = seqN(fix<nPED * 3 + 1>, fix<6>, fix<nPEM>);
        const auto &seqN_ColN3 = seqN(fix<nPED * 3 + 3>, fix<6>, fix<nPEM>);
        const auto &seqN_Col40 = seqN(fix<0>, fix<6>, fix<4>);
        const auto &seqN_Col41 = seqN(fix<1>, fix<6>, fix<4>);
        const auto &seqN_Col42 = seqN(fix<2>, fix<6>, fix<4>);
        const auto &seqN_Col43 = seqN(fix<3>, fix<6>, fix<4>);
        
        // subtract
        for (int isls = 0; isls < nsls(); isls++) {
            const eigen::RMatX46 &memVar4 = mMemVar_CD[isls];
            stressN(seqN_Row, seqN_ColN0) -= memVar4(seqN_Row, seqN_Col40);
            stressN(seqN_Row, seqN_ColN1) -= memVar4(seqN_Row, seqN_Col41);
            stressN(seqN_Row, seqN_ColN2) -= memVar4(seqN_Row, seqN_Col42);
            stressN(seqN_Row, seqN_ColN3) -= memVar4(seqN_Row, seqN_Col43);
        }
    }
    
    
    /////////////////// strain => strain4 ///////////////////
    // CASE == _1D_FR
    static void formStrain4(const eigen::vec_ar6_CMatPP_RM &strainN,
                            eigen::vec_ar6_CMat22_RM &strain4, int nx) {
        for (int idim = 0; idim < 6; idim++) {
            strain4[nx][idim](0, 0) = strainN[nx][idim](1, 1);
            strain4[nx][idim](0, 1) = strainN[nx][idim](1, 3);
            strain4[nx][idim](1, 0) = strainN[nx][idim](3, 1);
            strain4[nx][idim](1, 1) = strainN[nx][idim](3, 3);
        }
    }
    // CASE != _1D_FR
    static void formStrain4(const eigen::RMatXN6 &strainN,
                            eigen::RMatX46 &strain4, int nx) {
        using spectral::nPEM;
        using spectral::nPED;
        
        // seqN
        using Eigen::seqN;
        using Eigen::fix;
        const auto &seqN_Row = seqN(fix<0>, nx);
        const auto &seqN_ColN0 = seqN(fix<nPED * 1 + 1>, fix<6>, fix<nPEM>);
        const auto &seqN_ColN1 = seqN(fix<nPED * 1 + 3>, fix<6>, fix<nPEM>);
        const auto &seqN_ColN2 = seqN(fix<nPED * 3 + 1>, fix<6>, fix<nPEM>);
        const auto &seqN_ColN3 = seqN(fix<nPED * 3 + 3>, fix<6>, fix<nPEM>);
        const auto &seqN_Col40 = seqN(fix<0>, fix<6>, fix<4>);
        const auto &seqN_Col41 = seqN(fix<1>, fix<6>, fix<4>);
        const auto &seqN_Col42 = seqN(fix<2>, fix<6>, fix<4>);
        const auto &seqN_Col43 = seqN(fix<3>, fix<6>, fix<4>);
        
        // extract
        strain4(seqN_Row, seqN_Col40) = strainN(seqN_Row, seqN_ColN0);
        strain4(seqN_Row, seqN_Col41) = strainN(seqN_Row, seqN_ColN1);
        strain4(seqN_Row, seqN_Col42) = strainN(seqN_Row, seqN_ColN2);
        strain4(seqN_Row, seqN_Col43) = strainN(seqN_Row, seqN_ColN3);
    }
    
    
    /////////////////// strain4 => dStress ///////////////////
    template <CaseFA CASE, class FMat6> static void
    strainToStressCG4(const FMat6 &strain4, FMat6 &dStress, int nx,
                      const fa4::Property4 &dLambda,
                      const fa4::Property4 &dMu,
                      const fa4::Property4 &dMu2) {
        typedef fa4::FieldArithmetic4 FA;
        // use block 3 to store sii
        FA::F123xP<CASE>(nx, dStress, 3, strain4, 0, 1, 2, dLambda);
        // diagonal
        FA::R1_FxP<CASE>(nx, dStress, 0, 3, strain4, 0, dMu2);
        FA::R1_FxP<CASE>(nx, dStress, 1, 3, strain4, 1, dMu2);
        FA::R1_FxP<CASE>(nx, dStress, 2, 3, strain4, 2, dMu2);
        // off-diagonal
        FA::FxP<CASE>(nx, dStress, 3, strain4, 3, dMu);
        FA::FxP<CASE>(nx, dStress, 4, strain4, 4, dMu);
        FA::FxP<CASE>(nx, dStress, 5, strain4, 5, dMu);
    }
    
    
    /////////////////// apply ///////////////////
    // CASE == _1D_FR
    template <CaseFA CASE, class FMat6, class FCG6>
    void apply(const FMat6 &strain, FMat6 &stress, int nx,
               FCG6 &dStress, std::vector<FCG6> &memVar, FCG6 &strain4,
               const fa4::Property4 &dLambda,
               const fa4::Property4 &dMu,
               const fa4::Property4 &dMu2) const {
        // stress -= memory
        subtractMemFromStress(stress, nx);
        
        // update memory variables, phase 1
        Attenuation::updateMemAlphaBeta<CASE>(dStress, memVar, nx);
        
        // form strain CG4
        formStrain4(strain, strain4, nx);
        
        // update memory stress
        strainToStressCG4<CASE>(strain4, dStress, nx,
                                dLambda, dMu, dMu2);
        
        // update memory variables, phase 2
        Attenuation::updateMemGamma<CASE>(dStress, memVar, nx);
    }
};

#endif /* AttenuationCG4_hpp */
