//
//  Domain.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/1/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  computational domain

#include "Domain.hpp"
// mesh
#include "SolidPoint.hpp"
#include "FluidPoint.hpp"
#include "SolidElement.hpp"
#include "FluidElement.hpp"
// boundary
#include "SolidFluidBoundary.hpp"
#include "AbsorbingBoundary.hpp"
#include "AxialBoundary.hpp"
#include "FluidSurfaceBoundary.hpp"
// mpi
#include "Messaging.hpp"
// source
#include "WavefieldInjection.hpp"
#include "ElementSource.hpp"
// output
#include "StationGroup.hpp"
// verbose
#include "io.hpp"
#include "mpi.hpp"
#include "bstring.hpp"
#include "vector_tools.hpp"

// constructor
Domain::Domain() {
    // set all single unique_ptr's to nullptr
    // cannot do this in .hpp as it requires #include
    mSolidFluidBoundary = nullptr;
    mAxialBoundary = nullptr;
    mAbsorbingBoundary = nullptr;
    mFluidSurfaceBoundary = nullptr;
    mMessaging = nullptr;
    mWavefieldInjection = nullptr;
}

////////////////////// domain construction //////////////////////
// add a solid point
void Domain::addSolidPoint(const std::shared_ptr<SolidPoint> &point) {
    point->setDomainTag((int)mSolidPoints.size());
    mSolidPoints.push_back(point);
}

// add a fluid point
void Domain::addFluidPoint(const std::shared_ptr<FluidPoint> &point) {
    point->setDomainTag((int)mFluidPoints.size());
    mFluidPoints.push_back(point);
}

// add a solid element
void Domain::addSolidElement(const std::shared_ptr<SolidElement> &element) {
    element->setDomainTag((int)mSolidElements.size());
    mSolidElements.push_back(element);
}

// add a fluid element
void Domain::addFluidElement(const std::shared_ptr<FluidElement> &element) {
    element->setDomainTag((int)mFluidElements.size());
    mFluidElements.push_back(element);
}

// replace a solid element
void Domain::replaceSolidElement(const std::shared_ptr<SolidElement> &element) {
    // check domain tag
    int domainTag = element->getDomainTag();
    if (domainTag == -1) {
        throw std::runtime_error("Domain::replaceSolidElement || "
                                 "Domain tag must be set before replacement.");
    }
    // check quad tag
    if (element->getQuadTag() != mSolidElements[domainTag]->getQuadTag()) {
        throw std::runtime_error("Domain::replaceSolidElement || "
                                 "The new and the old elements must have "
                                 "the same Quad tag.");
    }
    //  check old reference
    if (mSolidElements[domainTag].use_count() > 1) {
        throw std::runtime_error("Domain::replaceSolidElement || "
                                 "The old element has been referred outside "
                                 "this domain and is thus irreplaceable.");
    }
    // after replacement, the old will be deleted
    mSolidElements[element->getDomainTag()] = element;
}

// replace a fluid element
void Domain::replaceFluidElement(const std::shared_ptr<FluidElement> &element) {
    // check domain tag
    int domainTag = element->getDomainTag();
    if (domainTag == -1) {
        throw std::runtime_error("Domain::replaceFluidElement || "
                                 "Domain tag must be set before replacement.");
    }
    // check quad tag
    if (element->getQuadTag() != mFluidElements[domainTag]->getQuadTag()) {
        throw std::runtime_error("Domain::replaceFluidElement || "
                                 "The new and the old elements must have "
                                 "the same Quad tag.");
    }
    //  check old reference
    if (mFluidElements[domainTag].use_count() > 1) {
        throw std::runtime_error("Domain::replaceFluidElement || "
                                 "The old element has been referred outside "
                                 "this domain and is thus irreplaceable.");
    }
    // after replacement, the old will be deleted
    mFluidElements[element->getDomainTag()] = element;
}

// set solid-fluid boundary
void Domain::
setSolidFluidBoundary(std::unique_ptr<const SolidFluidBoundary> &sfb) {
    mSolidFluidBoundary = std::move(sfb);
}

// set absorbing boundary
void Domain::
setAbsorbingBoundary(std::unique_ptr<const AbsorbingBoundary> &abc) {
    mAbsorbingBoundary = std::move(abc);
}

// set axial boundary
void Domain::
setAxialBoundary(std::unique_ptr<const AxialBoundary> &axb) {
    mAxialBoundary = std::move(axb);
}

// set fluid surface boundary
void Domain::
setFluidSurfaceBoundary(std::unique_ptr<const FluidSurfaceBoundary> &fsb) {
    mFluidSurfaceBoundary = std::move(fsb);
}

// set mpi messaging
void Domain::setMessaging(std::unique_ptr<Messaging> &msg) {
    mMessaging = std::move(msg);
}

// verbose
std::string Domain::verbose() const {
    //////////////////////////// count info ////////////////////////////
    // point
    std::map<std::string, int> typeCountPoint;
    for (const auto &point: mSolidPoints) {
        if (!mMessaging->pointInSmallerRank(point)) {
            vector_tools::aggregate(typeCountPoint, point->typeInfo(), 1);
        }
    }
    for (const auto &point: mFluidPoints) {
        if (!mMessaging->pointInSmallerRank(point)) {
            vector_tools::aggregate(typeCountPoint, point->typeInfo(), 1);
        }
    }
    mpi::aggregate(typeCountPoint, 0);
    
    // element
    std::map<std::string, int> typeCountElement;
    for (const auto &elem: mSolidElements) {
        vector_tools::aggregate(typeCountElement, elem->typeInfo(), 1);
    }
    for (const auto &elem: mFluidElements) {
        vector_tools::aggregate(typeCountElement, elem->typeInfo(), 1);
    }
    mpi::aggregate(typeCountElement, 0);
    
    // solid-fluid boundary
    std::map<std::string, int> typeCountSFB =
    mSolidFluidBoundary->countInfo(*mMessaging);
    mpi::aggregate(typeCountSFB, 0);
    
    // absorbing boundary
    std::map<std::string, int> typeCountABB =
    mAbsorbingBoundary->countInfo(*mMessaging);
    mpi::aggregate(typeCountABB, 0);
    
    // axial boundary
    int axSolid, axFluid;
    mAxialBoundary->countInfo(*mMessaging, axSolid, axFluid);
    axSolid = mpi::sum(axSolid);
    axFluid = mpi::sum(axFluid);
    
    // fluid surface boundary
    int fluidSurf = mFluidSurfaceBoundary->countInfo(*mMessaging);
    fluidSurf = mpi::sum(fluidSurf);
    
    // mpi
    int nRankComm = mMessaging->getNumRankComm();
    nRankComm = mpi::sum(nRankComm);
    nRankComm /= 2; // one comm counted by two ranks
    
    
    //////////////////////////// max key size ////////////////////////////
    int mkl = (int)std::string("rank-to-rank communication").size();
    mkl = std::max(mkl, vector_tools::maxKeyLength(typeCountPoint));
    mkl = std::max(mkl, vector_tools::maxKeyLength(typeCountElement));
    mkl = std::max(mkl, vector_tools::maxKeyLength(typeCountSFB));
    mkl = std::max(mkl, vector_tools::maxKeyLength(typeCountABB));
    
    
    //////////////////////////// write box ////////////////////////////
    std::stringstream ss;
    ss << bstring::boxTitle("Computational Domain");
    // point
    int total = vector_tools::sumValues(typeCountPoint);
    ss << bstring::boxSubTitle(0, "GLL points");
    ss << bstring::boxEquals(2, mkl, typeCountPoint);
    ss << bstring::boxEquals(2, mkl + 1, "Σ", total);
    // element
    total = vector_tools::sumValues(typeCountElement);
    ss << bstring::boxSubTitle(0, "Spectral elements");
    ss << bstring::boxEquals(2, mkl, typeCountElement);
    ss << bstring::boxEquals(2, mkl + 1, "Σ", total);
    // solid-fluid boundary
    total = vector_tools::sumValues(typeCountSFB);
    if (total > 0) {
        ss << bstring::boxSubTitle(0, "Solid-fluid boundary");
        ss << bstring::boxEquals(2, mkl, typeCountSFB);
        ss << bstring::boxEquals(2, mkl + 1, "Σ", total);
    }
    // absorbing boundary
    total = vector_tools::sumValues(typeCountABB);
    if (total > 0) {
        ss << bstring::boxSubTitle(0, "Absorbing boundary");
        ss << bstring::boxEquals(2, mkl, typeCountABB);
        ss << bstring::boxEquals(2, mkl + 1, "Σ", total);
    }
    // axial boundary
    ss << bstring::boxSubTitle(0, "Axial boundary");
    ss << bstring::boxEquals(2, mkl, "on solid domain", axSolid);
    ss << bstring::boxEquals(2, mkl, "on fluid domain", axFluid);
    ss << bstring::boxEquals(2, mkl + 1, "Σ ", axSolid + axFluid);
    // fluid surface boundary
    if (fluidSurf > 0) {
        ss << bstring::boxSubTitle(0, "Fluid surface boundary");
        ss << bstring::boxEquals(2, mkl + 1, "Σ", fluidSurf);
    }
    // mpi
    ss << bstring::boxSubTitle(0, "Domain decomposition");
    ss << bstring::boxEquals(2, mkl, "sub-domain (nproc)", mpi::nproc());
    ss << bstring::boxEquals(2, mkl, "rank-to-rank communication", nRankComm);
    // end
    ss << bstring::boxBaseline() << "\n\n";
    return ss.str();
}

// MESHING STAGE //
///////////////////
// SOLVING STAGE //

// set wavefield injection
void Domain::
setWavefieldInjection(std::unique_ptr<WavefieldInjection> &wj) {
    wj->setInDomain(*this);
    mWavefieldInjection = std::move(wj);
}

void Domain::
addElementSource(std::unique_ptr<const ElementSource> &esrc) {
    mElementSources.push_back(std::move(esrc));
}

// add station group in solid
void Domain::
addStationGroupInSolid(std::unique_ptr<StationGroupInSolid> &stgrp) {
    mStationGroupInSolids.push_back(std::move(stgrp));
}

// add station group in fluid
void Domain::
addStationGroupInFluid(std::unique_ptr<StationGroupInFluid> &stgrp) {
    mStationGroupInFluids.push_back(std::move(stgrp));
}


////////////////////// timeloop //////////////////////
// apply sources
void Domain::applySources(int tstep, double time) const {
    // wavefield injection
    if (mWavefieldInjection) {
        mWavefieldInjection->apply(tstep, time);
    }
    
    // general sources
    for (const std::unique_ptr<const ElementSource> &src: mElementSources) {
        src->apply(tstep, time);
    }
}

// compute stiffness
void Domain::computeStiffness() const {
    // solid
    for (const std::shared_ptr<const SolidElement> &element: mSolidElements) {
        element->computeStiff();
    }
    
    // fluid
    for (const std::shared_ptr<const FluidElement> &element: mFluidElements) {
        element->computeStiff();
    }
}

// apply solid-fluid boundary
void Domain::applySolidFluidBoundary() const {
    mSolidFluidBoundary->apply();
}

// apply point-wise boundary conditons
void Domain::applyPointBoundaries() const {
    // the following order matters!
    mAbsorbingBoundary->apply();
    mAxialBoundary->apply();
    mFluidSurfaceBoundary->apply();
}

// stiff to accel on points
void Domain::computeStiffToAccel() const {
    // solid
    for (const std::shared_ptr<SolidPoint> &point: mSolidPoints) {
        point->computeStiffToAccel();
    }
    
    // fluid
    for (const std::shared_ptr<FluidPoint> &point: mFluidPoints) {
        point->computeStiffToAccel();
    }
}

// mpi phase 1: commGatherSendRecv
void Domain::mpiGatherSendRecv() const {
    mMessaging->commGatherSendRecv();
}

// mpi phase 2: commWaitScatter
void Domain::mpiWaitScatter() const {
    mMessaging->commWaitScatter();
}

// initialize stations
void Domain::initializeStations() const {
    // solid
    for (const std::unique_ptr<StationGroupInSolid> &stgrp:
         mStationGroupInSolids) {
        stgrp->initialize();
    }
    
    // fluid
    for (const std::unique_ptr<StationGroupInFluid> &stgrp:
         mStationGroupInFluids) {
        stgrp->initialize();
    }
}

// record stations
void Domain::recordStations(int tstep, double time) const {
    // solid
    for (const std::unique_ptr<StationGroupInSolid> &stgrp:
         mStationGroupInSolids) {
        stgrp->record(tstep, time);
    }
    
    // fluid
    for (const std::unique_ptr<StationGroupInFluid> &stgrp:
         mStationGroupInFluids) {
        stgrp->record(tstep, time);
    }
}

// dump stations
void Domain::dumpStations() const {
    // solid
    for (const std::unique_ptr<StationGroupInSolid> &stgrp:
         mStationGroupInSolids) {
        stgrp->dumpToFile();
    }
    
    // fluid
    for (const std::unique_ptr<StationGroupInFluid> &stgrp:
         mStationGroupInFluids) {
        stgrp->dumpToFile();
    }
}

// finalize stations
void Domain::finalizeStations() const {
    // solid
    for (const std::unique_ptr<StationGroupInSolid> &stgrp:
         mStationGroupInSolids) {
        stgrp->finalize();
    }
    
    // fluid
    for (const std::unique_ptr<StationGroupInFluid> &stgrp:
         mStationGroupInFluids) {
        stgrp->finalize();
    }
}

// check stability
void Domain::checkStability(int tstep, double t, double dt) const {
    std::shared_ptr<Point> unstablePoint = nullptr;
    // solid
    for (const std::shared_ptr<SolidPoint> &point: mSolidPoints) {
        if (!point->stable()) {
            unstablePoint = point;
            break;
        }
    }
    // fluid
    if (!unstablePoint) {
        for (const std::shared_ptr<FluidPoint> &point: mFluidPoints) {
            if (!point->stable()) {
                unstablePoint = point;
                break;
            }
        }
    }
    // abort if unstable
    if (unstablePoint) {
        const eigen::DRow2 &sz = unstablePoint->getCoords() / 1e3;
        double r = sz.norm();
        double theta = (r < numerical::dEpsilon) ? 0. : acos(sz(1) / r);
        std::stringstream ss;
        ss << io::endl << bstring::boxTitle("AxiSEM3D ABORTED",
                                            bstring::fmt::gBoxWidth,
                                            bstring::fmt::gExceptionFill);
        ss << "SIMULATION BLEW UP!" << io::endl;
        ss << "Where the instability occured:" << io::endl;
        ss << bstring::boxEquals(2, 12, "(s, z) / km",
                                 bstring::range(sz(0), sz(1), '(', ')'));
        ss << bstring::boxEquals(2, 12, "radius / km", r);
        ss << bstring::boxEquals(2, 12, "theta / deg",
                                 theta / numerical::dDegree);
        ss << "When the instability occured:" << io::endl;
        ss << bstring::boxEquals(2, 12, "current time", t);
        ss << bstring::boxEquals(2, 12, "current step", tstep);
        ss << bstring::boxEquals(2, 12, "ΔT in use", dt);
        ss << bstring::boxBaseline(bstring::fmt::gBoxWidth,
                                   bstring::fmt::gExceptionFill) << io::endl;
        io::cout.setCoutWorldRank(mpi::rank());
        io::cout << ss.str();
        throw std::runtime_error("Domain::checkStability || "
                                 "Simulation blew up.");
    }
}


////////////////////////////////////////
//////////////// static ////////////////
////////////////////////////////////////

// stream-based creator
std::shared_ptr<Domain> Domain::
createFromStream(sbs::ifstream &ifs, const TimeScheme &timeScheme) {
    // create empty domain
    std::shared_ptr<Domain> domain = std::make_shared<Domain>();
    
    // solid points
    int nsp = ifs.get<int>();
    for (int ipnt = 0; ipnt < nsp; ipnt++) {
        ifs.get<std::string>(); // read class name "SolidPoint"
        domain->addSolidPoint(std::make_shared<SolidPoint>(ifs, timeScheme));
    }
    
    // fluid points
    int nfp = ifs.get<int>();
    for (int ipnt = 0; ipnt < nfp; ipnt++) {
        ifs.get<std::string>(); // read class name "FluidPoint"
        domain->addFluidPoint(std::make_shared<FluidPoint>(ifs, timeScheme));
    }
    
    // solid elements
    int nse = ifs.get<int>();
    for (int iele = 0; iele < nse; iele++) {
        ifs.get<std::string>(); // read class name "SolidElement"
        domain->addSolidElement(std::make_shared<SolidElement>(ifs, *domain));
    }
    
    // fluid elements
    int nfe = ifs.get<int>();
    for (int iele = 0; iele < nfe; iele++) {
        ifs.get<std::string>(); // read class name "FluidElement"
        domain->addFluidElement(std::make_shared<FluidElement>(ifs, *domain));
    }
    
    // solid-fluid boundary
    auto sfb = SolidFluidBoundary::createFromStream(ifs, *domain);
    domain->setSolidFluidBoundary(sfb);
    
    // absorbing boundary
    auto abb = AbsorbingBoundary::createFromStream(ifs, *domain);
    domain->setAbsorbingBoundary(abb);
    
    // axial boundary
    auto axb = AxialBoundary::createFromStream(ifs, *domain);
    domain->setAxialBoundary(axb);
    
    // fluid surface boundary
    auto fsb = FluidSurfaceBoundary::createFromStream(ifs, *domain);
    domain->setFluidSurfaceBoundary(fsb);
    
    // messaging
    auto msg = Messaging::createFromStream(ifs, *domain);
    domain->setMessaging(msg);
    
    // finish meshing stage
    return domain;
}
