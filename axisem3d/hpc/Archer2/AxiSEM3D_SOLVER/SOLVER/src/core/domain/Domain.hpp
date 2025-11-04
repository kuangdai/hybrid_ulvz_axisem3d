//
//  Domain.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/1/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  computational domain

#ifndef Domain_hpp
#define Domain_hpp

#include <vector>
#include <memory>

// mesh
class SolidPoint;
class FluidPoint;
class SolidElement;
class FluidElement;

// boundary
class SolidFluidBoundary;
class AbsorbingBoundary;
class AxialBoundary;
class FluidSurfaceBoundary;

// mpi
class Messaging;

// source
class WavefieldInjection;
class ElementSource;

// output
class StationSolid;
class StationFluid;
template <class StationT>
class StationGroup;
typedef StationGroup<StationSolid> StationGroupInSolid;
typedef StationGroup<StationFluid> StationGroupInFluid;

// stream
#include "sbstream.hpp"
class TimeScheme;

class Domain {
public:
    // constructor
    Domain();
    
    ////////////////////// domain construction //////////////////////
    // add a solid point
    void addSolidPoint(const std::shared_ptr<SolidPoint> &point);
    
    // add a fluid point
    void addFluidPoint(const std::shared_ptr<FluidPoint> &point);
    
    // add a solid element
    void addSolidElement(const std::shared_ptr<SolidElement> &element);
    
    // add a fluid element
    void addFluidElement(const std::shared_ptr<FluidElement> &element);
    
    // replace a solid element
    void replaceSolidElement(const std::shared_ptr<SolidElement> &element);
    
    // replace a fluid element
    void replaceFluidElement(const std::shared_ptr<FluidElement> &element);
    
    // set solid-fluid boundary
    void setSolidFluidBoundary(std::unique_ptr<const SolidFluidBoundary> &sfb);
    
    // set absorbing boundary
    void setAbsorbingBoundary(std::unique_ptr<const AbsorbingBoundary> &abc);
    
    // set axial boundary
    void setAxialBoundary(std::unique_ptr<const AxialBoundary> &axb);
    
    // set fluid surface boundary
    void
    setFluidSurfaceBoundary(std::unique_ptr<const FluidSurfaceBoundary> &fsb);
    
    // set mpi messaging
    void setMessaging(std::unique_ptr<Messaging> &msg);
    
    // verbose
    std::string verbose() const;
    
    // MESHING STAGE //
    ///////////////////
    // SOLVING STAGE //
    
    // set wavefield injection
    void setWavefieldInjection(std::unique_ptr<WavefieldInjection> &wj);
    
    // add a source
    void addElementSource(std::unique_ptr<const ElementSource> &esrc);
    
    // add station group in solid
    void addStationGroupInSolid(std::unique_ptr<StationGroupInSolid> &stgrp);
    
    // add station group in fluid
    void addStationGroupInFluid(std::unique_ptr<StationGroupInFluid> &stgrp);
    
    
    ////////////////////// get domain //////////////////////
    // get a solid point
    inline const std::shared_ptr<SolidPoint> &
    getSolidPoint(int domainTag) const {
        return mSolidPoints[domainTag];
    }
    
    // get a fluid point
    inline const std::shared_ptr<FluidPoint> &
    getFluidPoint(int domainTag) const {
        return mFluidPoints[domainTag];
    }
    
    // get solid points
    inline const std::vector<std::shared_ptr<SolidPoint>> &
    getSolidPoints() const {
        return mSolidPoints;
    }
    
    // get fluid points
    inline const std::vector<std::shared_ptr<FluidPoint>> &
    getFluidPoints() const {
        return mFluidPoints;
    }
    
    // get a solid element
    inline const std::shared_ptr<const SolidElement> &
    getSolidElement(int domainTag) const {
        return mSolidElements[domainTag];
    }
    
    // get a fluid element
    inline const std::shared_ptr<const FluidElement> &
    getFluidElement(int domainTag) const {
        return mFluidElements[domainTag];
    }
    
    // get solid elements
    inline const std::vector<std::shared_ptr<const SolidElement>> &
    getSolidElements() const {
        return mSolidElements;
    }
    
    // get fluid elements
    inline const std::vector<std::shared_ptr<const FluidElement>> &
    getFluidElements() const {
        return mFluidElements;
    }
    
    
    ////////////////////// timeloop //////////////////////
    // apply sources
    void applySources(int tstep, double time) const;
    
    // compute stiffness
    void computeStiffness() const;
    
    // apply solid-fluid boundary
    void applySolidFluidBoundary() const;
    
    // apply point-wise boundary conditons
    void applyPointBoundaries() const;
    
    // stiff to accel on points
    void computeStiffToAccel() const;
    
    // mpi phase 1: commGatherSendRecv
    void mpiGatherSendRecv() const;
    
    // mpi phase 2: commWaitScatter
    void mpiWaitScatter() const;
    
    // initialize stations
    void initializeStations() const;
    
    // record stations
    void recordStations(int tstep, double time) const;
    
    // dump stations
    void dumpStations() const;
    
    // finalize stations
    void finalizeStations() const;
    
    // check stability
    void checkStability(int tstep, double t, double dt) const;
    
    
private:
    ////////////////////// spectral elements //////////////////////
    // points
    std::vector<std::shared_ptr<SolidPoint>> mSolidPoints;
    std::vector<std::shared_ptr<FluidPoint>> mFluidPoints;
    
    // elements
    std::vector<std::shared_ptr<const SolidElement>> mSolidElements;
    std::vector<std::shared_ptr<const FluidElement>> mFluidElements;
    
    
    ////////////////////// boundaries //////////////////////
    // solid-fluid boundary
    std::unique_ptr<const SolidFluidBoundary> mSolidFluidBoundary;
    
    // axial boundary
    std::unique_ptr<const AxialBoundary> mAxialBoundary;
    
    // absorbing boundary
    std::unique_ptr<const AbsorbingBoundary> mAbsorbingBoundary;
    
    // fluid surface
    std::unique_ptr<const FluidSurfaceBoundary> mFluidSurfaceBoundary;
    
    
    ////////////////////// mpi //////////////////////
    std::unique_ptr<Messaging> mMessaging;
    
    // MESHING STAGE //
    ///////////////////
    // SOLVING STAGE //
    
    ////////////////////// source //////////////////////
    // injection
    std::unique_ptr<const WavefieldInjection> mWavefieldInjection;
    
    // sources
    std::vector<std::unique_ptr<const ElementSource>> mElementSources;
    
    
    ////////////////////// output //////////////////////
    // stations in solid
    std::vector<std::unique_ptr<StationGroupInSolid>> mStationGroupInSolids;
    
    // stations in fluid
    std::vector<std::unique_ptr<StationGroupInFluid>> mStationGroupInFluids;
    
    
    ////////////////////// learning //////////////////////
    // TODO
    
    
    ////////////////////////////////////////
    //////////////// static ////////////////
    ////////////////////////////////////////
    
public:
    // stream-based creator
    static std::shared_ptr<Domain>
    createFromStream(sbs::ifstream &ifs, const TimeScheme &timeScheme);
};

#endif /* Domain_hpp */
