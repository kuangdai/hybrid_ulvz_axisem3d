//
//  SolidElement.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/2/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  solid spectral element

#include "SolidElement.hpp"
#include "SolidPoint.hpp"
#include "fftw.hpp"
// output
#include "mapPPvsN.hpp"
// stream
#include "Domain.hpp"
// typeInfo
#include "bstring.hpp"

using spectral::nPEM;

// parameter constructor
SolidElement::SolidElement(int quadTag, std::unique_ptr<const GradQuad> &grad,
                           std::unique_ptr<const PRT> &prt,
                           std::unique_ptr<const Elastic> &elastic,
                           const std::array<std::shared_ptr<SolidPoint>,
                           spectral::nPEM> &points):
Element(quadTag, grad, prt), mElastic(elastic.release()),
mInFourier((mPRT ? mPRT->is1D() : true) && mElastic->is1D()),
mPoints(points) {
    // construct derived
    constructDerived();
}

// copy constructor
SolidElement::SolidElement(const SolidElement &other):
Element(other), mElastic(other.mElastic->clone()),
mInFourier((mPRT ? mPRT->is1D() : true) && mElastic->is1D()),
mPoints(other.mPoints) {
    // construct derived
    constructDerived();
}

// stream constructor
SolidElement::SolidElement(sbs::ifstream &ifs, const Domain &domain):
Element(ifs), mElastic(Elastic::createFromStream(ifs)),
mInFourier((mPRT ? mPRT->is1D() : true) && mElastic->is1D()),
mPoints() {
    // points
    const std::vector<int> &ptags = ifs.get<std::vector<int>>();
    for (int ipnt = 0; ipnt < nPEM; ipnt++) {
        mPoints[ipnt] = domain.getSolidPoint(ptags[ipnt]);
    }
    
    // construct derived
    constructDerived();
}

// construct derived
void SolidElement::constructDerived() {
    // point set
    pointSet(mInFourier);
    
    // check compatibility
    mElastic->checkCompatibility(mNr, mInFourier);
    
    // crd transform
    if (mElastic->inRTZ()) {
        setToRTZ_ByMaterialOrPRT();
    }
    
    // workspace
    if (sStrainSpherical_CD.rows() < mNr) {
        expandWorkspace(mNr);
    }
    
    // report request to FFTW
    if (!mInFourier) {
        if (mPRT) {
            fftw::gFFT_N9.addNR(mNr);
        } else {
            fftw::gFFT_N6.addNR(mNr);
        }
    }
}

// type info
std::string SolidElement::typeInfo() const {
    std::string info = "SolidElement";
    if (mPRT) {
        if (mPRT->is1D()) {
            info = info + "$PRT1D";
        } else {
            info = info + "$PRT3D";
        }
    }
    if (mElastic->is1D()) {
        info = info + "$" + bstring::typeName(*mElastic) + "1D";
    } else {
        info = info + "$" + bstring::typeName(*mElastic) + "3D";
    }
    return info;
}


/////////////////////////// point ///////////////////////////
// get point
Point &SolidElement::getPoint(int ipnt) const {
    return *(mPoints[ipnt]);
}


/////////////////////////// time loop ///////////////////////////
// collect displacement from points
void SolidElement::
collectDisplFromPoints(eigen::vec_ar3_CMatPP_RM &displElem) const {
    int nu_1 = mNu + 1;
    for (int ipol = 0; ipol < spectral::nPED; ipol++) {
        for (int jpol = 0; jpol < spectral::nPED; jpol++) {
            mPoints[ipol * spectral::nPED + jpol]->
            scatterDisplToElement(displElem, nu_1, ipol, jpol);
        }
    }
}

// displacement to stiffness
void SolidElement::
displToStiff(const eigen::vec_ar3_CMatPP_RM &displElem,
             eigen::vec_ar3_CMatPP_RM &stiffElem) const {
    int nu_1 = mNu + 1;
    if (mPRT) {
        /////////////// with PRT, strain in 3 * 3 ///////////////
        // disp to strain
        mGradQuad->computeGrad9(displElem, sStrainSpherical_FR, nu_1);
        mTransform->transformSPZ_RTZ9(sStrainSpherical_FR, nu_1);
        
        // strain to stress
        if (mInFourier) {
            mPRT->sphericalToUndulated6_FR(sStrainSpherical_FR,
                                           sStrainUndulated_FR, nu_1);
            mElastic->strainToStress_FR(sStrainUndulated_FR,
                                        sStressUndulated_FR, nu_1);
            // record stress if needed for output
            if (mStressBuffer->rows() > 0) {
                mapPPvsN::PP2N(sStressUndulated_FR, *mStressBuffer, nu_1);
            }
            mPRT->undulatedToSpherical6_FR(sStressUndulated_FR,
                                           sStressSpherical_FR, nu_1);
        } else {
            fftw::gFFT_N9.computeC2R(sStrainSpherical_FR,
                                     sStrainSpherical_CD, mNr);
            mPRT->sphericalToUndulated6_CD(sStrainSpherical_CD,
                                           sStrainUndulated_CD, mNr);
            mElastic->strainToStress_CD(sStrainUndulated_CD,
                                        sStressUndulated_CD, mNr);
            // record stress if needed for output
            if (mStressBuffer->rows() > 0) {
                /////////// additional support required from gFFT_N6 ///////////
                fftw::gFFT_N6.computeR2C(sStressUndulated_CD,
                                         sStressUndulated_FR, mNr);
                mapPPvsN::PP2N(sStressUndulated_FR, *mStressBuffer, nu_1);
            }
            mPRT->undulatedToSpherical6_CD(sStressUndulated_CD,
                                           sStressSpherical_CD, mNr);
            fftw::gFFT_N9.computeR2C(sStressSpherical_CD,
                                     sStressSpherical_FR, mNr);
        }
        
        // stress to stiffness
        mTransform->transformRTZ_SPZ9(sStressSpherical_FR, nu_1);
        mGradQuad->computeQuad9(sStressSpherical_FR, stiffElem, nu_1);
    } else {
        /////////////// without PRT, strain in 6 * 1, Voigt ///////////////
        // disp to strain
        mGradQuad->computeGrad6(displElem, sStrainUndulated_FR, nu_1);
        if (mElastic->inRTZ()) {
            mTransform->transformSPZ_RTZ6(sStrainUndulated_FR, nu_1);
        }
        
        // strain to stress
        if (mInFourier) {
            mElastic->strainToStress_FR(sStrainUndulated_FR,
                                        sStressUndulated_FR, nu_1);
        } else {
            fftw::gFFT_N6.computeC2R(sStrainUndulated_FR,
                                     sStrainUndulated_CD, mNr);
            mElastic->strainToStress_CD(sStrainUndulated_CD,
                                        sStressUndulated_CD, mNr);
            fftw::gFFT_N6.computeR2C(sStressUndulated_CD,
                                     sStressUndulated_FR, mNr);
        }
        
        // record stress if needed for output
        // record before rotation to keep RTZ consistent with strain
        if (mStressBuffer->rows() > 0) {
            mapPPvsN::PP2N(sStressUndulated_FR, *mStressBuffer, nu_1);
        }
        
        // stress to stiffness
        if (mElastic->inRTZ()) {
            mTransform->transformRTZ_SPZ6(sStressUndulated_FR, nu_1);
        }
        mGradQuad->computeQuad6(sStressUndulated_FR, stiffElem, nu_1);
    }
}

// add stiffness to points
// allow a derived class to change stiffElem (no const)
void SolidElement::
addStiffToPoints(eigen::vec_ar3_CMatPP_RM &stiffElem) const {
    for (int ipol = 0; ipol < spectral::nPED; ipol++) {
        for (int jpol = 0; jpol < spectral::nPED; jpol++) {
            mPoints[ipol * spectral::nPED + jpol]->
            gatherStiffFromElement(stiffElem, ipol, jpol);
        }
    }
}

// compute stiffness term
void SolidElement::computeStiff() const {
    // collect displacement from points
    collectDisplFromPoints(sDisplSpherical_FR);
    
    // displacement to stiffness
    displToStiff(sDisplSpherical_FR, sStiffSpherical_FR);
    
    // add stiffness to points
    addStiffToPoints(sStiffSpherical_FR);
}


/////////////////////////// source ///////////////////////////
// prepare force source (force given in SPZ)
void SolidElement::prepareForceSource() const {
    // seems nothing
}

// add force source (force given in SPZ)
void SolidElement::addForceSource(const eigen::CMatXN3 &force,
                                  int nu_1_force) const {
    for (int ipnt = 0; ipnt < nPEM; ipnt++) {
        mPoints[ipnt]->addForceSource(force, nu_1_force, ipnt);
    }
}

// prepare moment source (moment tensor given in SPZ)
void SolidElement::prepareMomentSource() const {
    // requires FFT N6 even with PRT
    if (mPRT) {
        if (!(mPRT->is1D())) {
            fftw::gFFT_N6.addNR(mNr);
        }
    }
}

// add moment source (moment tensor given in SPZ)
void SolidElement::addMomentSource(const eigen::CMatXN6 &moment,
                                   int nu_1_moment) const {
    // pad source with zeros if source has lower order than element
    // truncate source if source has higher order than element
    int nu_1_element = mNu + 1;
    int nu_1_coexist = std::min(nu_1_element, nu_1_moment);
    mapPPvsN::N2PP(moment, sStressUndulated_FR, nu_1_coexist);
    for (int alpha = nu_1_coexist; alpha < nu_1_element; alpha++) {
        for (int idim = 0; idim < 6; idim++) {
            sStressUndulated_FR[alpha][idim].setZero();
        }
    }
    
    // by default, source order does not change without 3D PRT
    int nu_1_source = nu_1_coexist;
    
    // stress to stiffness
    if (mPRT) {
        /////////////// with PRT, strain in 3 * 3 ///////////////
        // change to Voigt convention for strain before rotation
        for (int alpha = 0; alpha < nu_1_source; alpha++) {
            // dims 3, 4, 5 are shear components
            for (int idim = 3; idim < 6; idim++) {
                sStressUndulated_FR[alpha][idim] *= (numerical::Real)2.;
            }
        }
        mTransform->transformSPZ_RTZ6(sStressUndulated_FR, nu_1_source);
        for (int alpha = 0; alpha < nu_1_source; alpha++) {
            // dims 3, 4, 5 are shear components
            for (int idim = 3; idim < 6; idim++) {
                sStressUndulated_FR[alpha][idim] *= (numerical::Real).5;
            }
        }
        
        // PRT
        if (mPRT->is1D()) {
            mPRT->undulatedToSpherical6_NoIntegration_FR(sStressUndulated_FR,
                                                         sStressSpherical_FR,
                                                         nu_1_source);
        } else {
            fftw::gFFT_N6.computeC2R(sStressUndulated_FR,
                                     sStressUndulated_CD, mNr);
            mPRT->undulatedToSpherical6_NoIntegration_CD(sStressUndulated_CD,
                                                         sStressSpherical_CD,
                                                         mNr);
            fftw::gFFT_N9.computeR2C(sStressSpherical_CD,
                                     sStressSpherical_FR, mNr);
            // source order may increase with 3D PRT
            nu_1_source = nu_1_element;
        }
        
        // stress to stiffness
        mTransform->transformRTZ_SPZ9(sStressSpherical_FR, nu_1_source);
        mGradQuad->computeQuad9_NoIntegration(sStressSpherical_FR,
                                              sStiffSpherical_FR, nu_1_source);
    } else {
        /////////////// without PRT, strain in 6 * 1, Voigt ///////////////
        // stress to stiffness
        mGradQuad->computeQuad6_NoIntegration(sStressUndulated_FR,
                                              sStiffSpherical_FR,
                                              nu_1_source);
    }
    
    // add stiffness to points
    // NOTE: doing the following with addStiffToPoints(sStiffSpherical_FR)
    // is incorrect if the moment tensor is on the injection boundary
    for (int ipol = 0; ipol < spectral::nPED; ipol++) {
        for (int jpol = 0; jpol < spectral::nPED; jpol++) {
            mPoints[ipol * spectral::nPED + jpol]->
            gatherStiffFromElement(sStiffSpherical_FR, ipol, jpol);
        }
    }
}


/////////////////////////// wavefield output ///////////////////////////
// prepare wavefield output
void SolidElement::
prepareWavefieldOutput(const channel::solid::ChannelOptions &chops) const {
    // buffer
    if (chops.mNeedBufferS) {
        mStressBuffer->resize(mNu + 1, spectral::nPEM * 6);
    }
    
    // fft
    if (mPRT && !mInFourier && (chops.mNeedBufferE || chops.mNeedBufferS)) {
        fftw::gFFT_N6.addNR(mNr);
    }
}

// displ field
void SolidElement::getDisplField(eigen::CMatXN3 &displ) const {
    // collect displacement from points
    collectDisplFromPoints(sDisplSpherical_FR);
    
    // copy
    mapPPvsN::PP2N(sDisplSpherical_FR, displ, mNu + 1);
}

// nabla field
void SolidElement::getNablaField(eigen::CMatXN9 &nabla) const {
    // collect displacement from points
    collectDisplFromPoints(sDisplSpherical_FR);
    int nu_1 = mNu + 1;
    
    // displ to nabla, use sStressSpherical as temp storage
    eigen::vec_ar9_CMatPP_RM &sStrainUndulated9_FR = sStressSpherical_FR;
    eigen::RMatXN9 &sStrainUndulated9_CD = sStressSpherical_CD;
    if (mPRT) {
        mGradQuad->computeGrad9(sDisplSpherical_FR,
                                sStrainSpherical_FR, nu_1);
        mTransform->transformSPZ_RTZ9(sStrainSpherical_FR, nu_1);
        if (mInFourier) {
            mPRT->sphericalToUndulated9_FR(sStrainSpherical_FR,
                                           sStrainUndulated9_FR, nu_1);
        } else {
            fftw::gFFT_N9.computeC2R(sStrainSpherical_FR,
                                     sStrainSpherical_CD, mNr);
            mPRT->sphericalToUndulated9_CD(sStrainSpherical_CD,
                                           sStrainUndulated9_CD, nu_1);
            fftw::gFFT_N9.computeR2C(sStrainUndulated9_CD,
                                     sStrainUndulated9_FR, mNr);
        }
    } else {
        mGradQuad->computeGrad9(sDisplSpherical_FR,
                                sStrainUndulated9_FR, nu_1);
    }
    
    // copy
    mapPPvsN::PP2N(sStrainUndulated9_FR, nabla, nu_1);
}

// strain field
void SolidElement::getStrainField(eigen::CMatXN6 &strain) const {
    // collect displacement from points
    collectDisplFromPoints(sDisplSpherical_FR);
    int nu_1 = mNu + 1;
    
    // displ to strain
    if (mPRT) {
        mGradQuad->computeGrad9(sDisplSpherical_FR,
                                sStrainSpherical_FR, nu_1);
        mTransform->transformSPZ_RTZ9(sStrainSpherical_FR, nu_1);
        if (mInFourier) {
            mPRT->sphericalToUndulated6_FR(sStrainSpherical_FR,
                                           sStrainUndulated_FR, nu_1);
        } else {
            fftw::gFFT_N9.computeC2R(sStrainSpherical_FR,
                                     sStrainSpherical_CD, mNr);
            mPRT->sphericalToUndulated6_CD(sStrainSpherical_CD,
                                           sStrainUndulated_CD, mNr);
            /////////// additional support required from gFFT_N6 ///////////
            fftw::gFFT_N6.computeR2C(sStrainUndulated_CD,
                                     sStrainUndulated_FR, mNr);
        }
    } else {
        mGradQuad->computeGrad6(sDisplSpherical_FR,
                                sStrainUndulated_FR, nu_1);
    }
    
    // copy
    mapPPvsN::PP2N(sStrainUndulated_FR, strain, nu_1);
}

// curl field
void SolidElement::getCurlField(eigen::CMatXN3 &curl) const {
    // collect displacement from points
    collectDisplFromPoints(sDisplSpherical_FR);
    int nu_1 = mNu + 1;
    
    // displ to nabla, use sStressSpherical as temp storage
    eigen::vec_ar9_CMatPP_RM &sStrainUndulated9_FR = sStressSpherical_FR;
    eigen::RMatXN9 &sStrainUndulated9_CD = sStressSpherical_CD;
    if (mPRT) {
        mGradQuad->computeGrad9(sDisplSpherical_FR,
                                sStrainSpherical_FR, nu_1);
        mTransform->transformSPZ_RTZ9(sStrainSpherical_FR, nu_1);
        if (mInFourier) {
            mPRT->sphericalToUndulated9_FR(sStrainSpherical_FR,
                                           sStrainUndulated9_FR, nu_1);
        } else {
            fftw::gFFT_N9.computeC2R(sStrainSpherical_FR,
                                     sStrainSpherical_CD, mNr);
            mPRT->sphericalToUndulated9_CD(sStrainSpherical_CD,
                                           sStrainUndulated9_CD, nu_1);
            fftw::gFFT_N9.computeR2C(sStrainUndulated9_CD,
                                     sStrainUndulated9_FR, mNr);
        }
    } else {
        mGradQuad->computeGrad9(sDisplSpherical_FR,
                                sStrainUndulated9_FR, nu_1);
    }
    
    // nabla to curl, use sStiffSpherical_FR as temp storage
    eigen::vec_ar3_CMatPP_RM &sCurlUndulated_FR = sStiffSpherical_FR;
    for (int alpha = 0; alpha < nu_1; alpha++) {
        sCurlUndulated_FR[alpha][0] = (sStrainUndulated9_FR[alpha][7] -
                                       sStrainUndulated9_FR[alpha][5]);
        sCurlUndulated_FR[alpha][1] = (sStrainUndulated9_FR[alpha][2] -
                                       sStrainUndulated9_FR[alpha][6]);
        sCurlUndulated_FR[alpha][2] = (sStrainUndulated9_FR[alpha][3] -
                                       sStrainUndulated9_FR[alpha][1]);
    }
    
    // copy
    mapPPvsN::PP2N(sCurlUndulated_FR, curl, nu_1);
}
