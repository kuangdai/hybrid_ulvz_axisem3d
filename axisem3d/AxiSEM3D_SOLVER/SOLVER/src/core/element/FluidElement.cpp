//
//  FluidElement.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/2/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  fluid spectral element

#include "FluidElement.hpp"
#include "FluidPoint.hpp"
#include "fftw.hpp"
// output
#include "mapPPvsN.hpp"
// stream
#include "Domain.hpp"

using spectral::nPEM;

// parameter constructor
FluidElement::FluidElement(int quadTag, std::unique_ptr<const GradQuad> &grad,
                           std::unique_ptr<const PRT> &prt,
                           std::unique_ptr<const Acoustic> &acoustic,
                           const std::array<std::shared_ptr<FluidPoint>,
                           spectral::nPEM> &points):
Element(quadTag, grad, prt), mAcoustic(acoustic.release()),
mInFourier((mPRT ? mPRT->is1D() : true) && mAcoustic->is1D()),
mPoints(points) {
    // construct derived
    constructDerived();
}

// copy constructor
FluidElement::FluidElement(const FluidElement &other):
Element(other), mAcoustic(std::make_unique<Acoustic>(*(other.mAcoustic))),
mInFourier((mPRT ? mPRT->is1D() : true) && mAcoustic->is1D()),
mPoints(other.mPoints) {
    // construct derived
    constructDerived();
}

// stream constructor
FluidElement::FluidElement(sbs::ifstream &ifs, const Domain &domain):
Element(ifs), mAcoustic(std::make_unique<Acoustic>(ifs)),
mInFourier((mPRT ? mPRT->is1D() : true) && mAcoustic->is1D()),
mPoints() {
    // points
    const std::vector<int> &ptags = ifs.get<std::vector<int>>();
    for (int ipnt = 0; ipnt < nPEM; ipnt++) {
        mPoints[ipnt] = domain.getFluidPoint(ptags[ipnt]);
    }
    
    // construct derived
    constructDerived();
}

// construct derived
void FluidElement::constructDerived() {
    // point set
    pointSet(mInFourier);
    
    // check compatibility
    mAcoustic->checkCompatibility(mNr, mInFourier);
    
    // workspace
    if (sStrainSpherical_CD.rows() < mNr) {
        expandWorkspace(mNr);
    }
    
    // report request to FFTW
    if (!mInFourier) {
        fftw::gFFT_N3.addNR(mNr);
    }
}

// type info
std::string FluidElement::typeInfo() const {
    std::string info = "FluidElement";
    if (mPRT) {
        if (mPRT->is1D()) {
            info = info + "$PRT1D";
        } else {
            info = info + "$PRT3D";
        }
    }
    if (mAcoustic->is1D()) {
        info = info + "$Acoustic1D";
    } else {
        info = info + "$Acoustic3D";
    }
    return info;
}


/////////////////////////// point ///////////////////////////
// get point
Point &FluidElement::getPoint(int ipnt) const {
    return *(mPoints[ipnt]);
}


/////////////////////////// time loop ///////////////////////////
// collect displacement from points
void FluidElement::
collectDisplFromPoints(eigen::vec_ar1_CMatPP_RM &displElem) const {
    int nu_1 = mNu + 1;
    for (int ipol = 0; ipol < spectral::nPED; ipol++) {
        for (int jpol = 0; jpol < spectral::nPED; jpol++) {
            mPoints[ipol * spectral::nPED + jpol]->
            scatterDisplToElement(displElem, nu_1, ipol, jpol);
        }
    }
}

// displacement to stiffness
void FluidElement::
displToStiff(const eigen::vec_ar1_CMatPP_RM &displElem,
             eigen::vec_ar1_CMatPP_RM &stiffElem) const {
    int nu_1 = mNu + 1;
    
    // displ => strain
    mGradQuad->computeGrad3(displElem,
                            sStrainSpherical_FR, nu_1);
    
    // strain => stress
    if (mPRT) {
        mTransform->transformSPZ_RTZ3(sStrainSpherical_FR, nu_1);
        if (mInFourier) {
            mPRT->sphericalToUndulated3_FR(sStrainSpherical_FR,
                                           sStrainUndulated_FR, nu_1);
            mAcoustic->strainToStress_FR(sStrainUndulated_FR,
                                         sStressUndulated_FR, nu_1);
            mPRT->undulatedToSpherical3_FR(sStressUndulated_FR,
                                           sStressSpherical_FR, nu_1);
        } else {
            fftw::gFFT_N3.computeC2R(sStrainSpherical_FR,
                                     sStrainSpherical_CD, mNr);
            mPRT->sphericalToUndulated3_CD(sStrainSpherical_CD,
                                           sStrainUndulated_CD, mNr);
            mAcoustic->strainToStress_CD(sStrainUndulated_CD,
                                         sStressUndulated_CD, mNr);
            mPRT->undulatedToSpherical3_CD(sStressUndulated_CD,
                                           sStressSpherical_CD, mNr);
            fftw::gFFT_N3.computeR2C(sStressSpherical_CD,
                                     sStressSpherical_FR, mNr);
        }
        mTransform->transformRTZ_SPZ3(sStressSpherical_FR, nu_1);
    } else {
        if (mInFourier) {
            mAcoustic->strainToStress_FR(sStrainSpherical_FR,
                                         sStressSpherical_FR, nu_1);
        } else {
            fftw::gFFT_N3.computeC2R(sStrainSpherical_FR,
                                     sStrainSpherical_CD, mNr);
            mAcoustic->strainToStress_CD(sStrainSpherical_CD,
                                         sStressSpherical_CD, mNr);
            fftw::gFFT_N3.computeR2C(sStressSpherical_CD,
                                     sStressSpherical_FR, mNr);
        }
    }
    
    // stress => stiffness
    mGradQuad->computeQuad3(sStressSpherical_FR,
                            stiffElem, nu_1);
}

// add stiffness to points
// allow a derived class to change stiffElem (no const)
void FluidElement::
addStiffToPoints(eigen::vec_ar1_CMatPP_RM &stiffElem) const {
    for (int ipol = 0; ipol < spectral::nPED; ipol++) {
        for (int jpol = 0; jpol < spectral::nPED; jpol++) {
            mPoints[ipol * spectral::nPED + jpol]->
            gatherStiffFromElement(stiffElem, ipol, jpol);
        }
    }
}

/////////////////////////// time loop ///////////////////////////
// compute stiffness term
void FluidElement::computeStiff() const {
    // collect displacement from points
    collectDisplFromPoints(sDisplSpherical_FR);
    
    // displacement to stiffness
    displToStiff(sDisplSpherical_FR, sStiffSpherical_FR);
    
    // add stiffness to points
    addStiffToPoints(sStiffSpherical_FR);
}


/////////////////////////// source ///////////////////////////
// prepare pressure source
void FluidElement::preparePressureSource() const {
    for (int ipnt = 0; ipnt < nPEM; ipnt++) {
        mPoints[ipnt]->preparePressureSource();
    }
}

// add pressure source
void FluidElement::addPressureSource(const eigen::CMatXN &pressure,
                                     int nu_1_pressure) const {
    for (int ipnt = 0; ipnt < nPEM; ipnt++) {
        mPoints[ipnt]->addPressureSource(pressure, nu_1_pressure, ipnt);
    }
}


/////////////////////////// wavefield output ///////////////////////////
// prepare wavefield output
void FluidElement::
prepareWavefieldOutput(const channel::fluid::ChannelOptions &chops) const {
    if (chops.mNeedBufferP) {
        for (int ipnt = 0; ipnt < nPEM; ipnt++) {
            mPoints[ipnt]->preparePressureOutput();
        }
    }
    
    if (chops.mNeedBufferD) {
        for (int ipnt = 0; ipnt < nPEM; ipnt++) {
            mPoints[ipnt]->prepareDeltaOutput();
        }
    }
}

// chi field
void FluidElement::getChiField(eigen::CMatXN &chi) const {
    // collect displacement from points
    collectDisplFromPoints(sDisplSpherical_FR);
    
    // convert to flattened
    mapPPvsN::PP2N(sDisplSpherical_FR, chi, mNu + 1);
}

// displ field
void FluidElement::getDisplField(eigen::CMatXN3 &displ) const {
    // collect displacement from points
    collectDisplFromPoints(sDisplSpherical_FR);
    int nu_1 = mNu + 1;
    
    ///////////////////// compute stress /////////////////////
    // displ => strain
    mGradQuad->computeGrad3(sDisplSpherical_FR,
                            sStrainSpherical_FR, nu_1);
    
    // strain => stress
    if (mPRT) {
        mTransform->transformSPZ_RTZ3(sStrainSpherical_FR, nu_1);
        if (mInFourier) {
            mPRT->sphericalToUndulated3_FR(sStrainSpherical_FR,
                                           sStrainUndulated_FR, nu_1);
            mAcoustic->strainToStress_FR(sStrainUndulated_FR,
                                         sStressUndulated_FR, nu_1);
        } else {
            fftw::gFFT_N3.computeC2R(sStrainSpherical_FR,
                                     sStrainSpherical_CD, mNr);
            mPRT->sphericalToUndulated3_CD(sStrainSpherical_CD,
                                           sStrainUndulated_CD, mNr);
            mAcoustic->strainToStress_CD(sStrainUndulated_CD,
                                         sStressUndulated_CD, mNr);
            fftw::gFFT_N3.computeR2C(sStressUndulated_CD,
                                     sStressUndulated_FR, mNr);
        }
        // convert to flattened
        mapPPvsN::PP2N(sStressUndulated_FR, displ, nu_1);
    } else {
        if (mInFourier) {
            mAcoustic->strainToStress_FR(sStrainSpherical_FR,
                                         sStressSpherical_FR, nu_1);
        } else {
            fftw::gFFT_N3.computeC2R(sStrainSpherical_FR,
                                     sStrainSpherical_CD, mNr);
            mAcoustic->strainToStress_CD(sStrainSpherical_CD,
                                         sStressSpherical_CD, mNr);
            fftw::gFFT_N3.computeR2C(sStressSpherical_CD,
                                     sStressSpherical_FR, mNr);
        }
        // convert to flattened
        mapPPvsN::PP2N(sStressSpherical_FR, displ, nu_1);
    }
}

// pressure field
void FluidElement::getPressureField(eigen::CMatXN &pressure) const {
    int nu_1 = mNu + 1;
    for (int ipnt = 0; ipnt < nPEM; ipnt++) {
        mPoints[ipnt]->scatterPressureToElement(pressure, nu_1, ipnt);
    }
}

// delta field
void FluidElement::getDeltaField(eigen::CMatXN &delta) const {
    int nu_1 = mNu + 1;
    for (int ipnt = 0; ipnt < nPEM; ipnt++) {
        mPoints[ipnt]->scatterDeltaToElement(delta, nu_1, ipnt);
    }
}

