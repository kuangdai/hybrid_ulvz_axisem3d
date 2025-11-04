//
//  NewmarkTimeScheme.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/30/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  Newmark time scheme

#include "NewmarkTimeScheme.hpp"
#include "SolidPoint.hpp"
#include "FluidPoint.hpp"
#include "Domain.hpp"
#include "timer.hpp"

//////////////////// point ////////////////////
// create fields on a solid point
void NewmarkTimeScheme::createFields(SolidPoint &point) const {
    int nu_1 = point.getNu() + 1;
    SolidPoint::Fields &f = point.getFields();
    f.mStiff = eigen::CMatX3::Zero(nu_1, 3);
    f.mDispl = eigen::CMatX3::Zero(nu_1, 3);
    f.mVeloc = eigen::CMatX3::Zero(nu_1, 3);
    f.mAccel = eigen::CMatX3::Zero(nu_1, 3);
}

// create fields on a fluid point
void NewmarkTimeScheme::createFields(FluidPoint &point) const {
    int nu_1 = point.getNu() + 1;
    FluidPoint::Fields &f = point.getFields();
    f.mStiff = eigen::CColX::Zero(nu_1, 1);
    f.mDispl = eigen::CColX::Zero(nu_1, 1);
    f.mVeloc = eigen::CColX::Zero(nu_1, 1);
    f.mAccel = eigen::CColX::Zero(nu_1, 1);
}


//////////////////// timeloop ////////////////////
// solve
void NewmarkTimeScheme::solve() const {
    // loop timer
    std::vector<std::string> timerNames = {"TOTAL", "SOURCE",
        "STIFFNESS", "MASS_TERM", "TIME_MARCH", "SF_COUPLING",
        "BOUNDARIES", "MPI_COMM", "WAVE_OUTPUT"};
    std::map<std::string, SimpleTimer> timers;
    for (const std::string &name: timerNames) {
        timers.insert({name, SimpleTimer()});
    }
    
    // verbose
    verboseBegin("NEWMARK");
    
    // start timer
    timers.at("TOTAL").resume();
    
    // wavefield output: initialize
    timers.at("WAVE_OUTPUT").resume();
    mDomain->initializeStations();
    timers.at("WAVE_OUTPUT").pause();
    
    // simulation time
    double t = mT0;
    
    ////////////////////////// loop //////////////////////////
    for (int tstep = 0; tstep < mNumTimeSteps; tstep++) {
        // apply source
        timers.at("SOURCE").resume();
        mDomain->applySources(tstep, t);
        timers.at("SOURCE").pause();
        
        // compute stiffness
        timers.at("STIFFNESS").resume();
        mDomain->computeStiffness();
        timers.at("STIFFNESS").pause();
        
        // solid-fluid boundary
        // must be applied after source and stiffness
        timers.at("SF_COUPLING").resume();
        mDomain->applySolidFluidBoundary();
        timers.at("SF_COUPLING").pause();
        
        // assemble phase 1: gather + send + recv
        timers.at("MPI_COMM").resume();
        mDomain->mpiGatherSendRecv();
        timers.at("MPI_COMM").pause();
        
        // wavefield output: record and dump
        timers.at("WAVE_OUTPUT").resume();
        mDomain->recordStations(tstep, t);
        timers.at("WAVE_OUTPUT").pause();
        
        // stability
        if (tstep % mStabilityInterval == 0) {
            mDomain->checkStability(tstep + 1, t, mDt);
        }
        
        // verbose
        verboseIter(timers.at("TOTAL").elapsedTotal(), tstep + 1, t);
        
        // assemble phase 2: wait + scatter
        timers.at("MPI_COMM").resume();
        mDomain->mpiWaitScatter();
        timers.at("MPI_COMM").pause();
        
        // apply point-wise boundary conditions
        timers.at("BOUNDARIES").resume();
        mDomain->applyPointBoundaries();
        timers.at("BOUNDARIES").pause();
        
        // stiff => accel
        timers.at("MASS_TERM").resume();
        mDomain->computeStiffToAccel();
        timers.at("MASS_TERM").pause();
        
        // update fields
        timers.at("TIME_MARCH").resume();
        update(mDomain->getSolidPoints());
        update(mDomain->getFluidPoints());
        timers.at("TIME_MARCH").pause();
        
        // march
        t += mDt;
    } // end of time loop
    
    // wavefield output: dump remaining buffers and finalize
    timers.at("WAVE_OUTPUT").resume();
    mDomain->dumpStations();
    mDomain->finalizeStations();
    timers.at("WAVE_OUTPUT").pause();
    
    // total
    timers.at("TOTAL").pause();
    
    // verbose
    verboseEnd("NEWMARK", timers.at("TOTAL").elapsedTotal(), t);
    
    // report timers
    reportLoopTimers(timerNames, timers);
}

// update fields on a solid point
void NewmarkTimeScheme::
update(const std::vector<std::shared_ptr<SolidPoint>> &points) const {
    static const numerical::Real dt = mDt;
    static const numerical::Real half_dt = .5 * mDt;
    static const numerical::Real half_dt_dt = .5 * mDt * mDt;
    for (const std::shared_ptr<SolidPoint> &sp: points) {
        // update dt
        SolidPoint::Fields &f = sp->getFields();
        f.mVeloc += half_dt * (f.mAccel + f.mStiff);
        f.mAccel = f.mStiff;
        f.mDispl += dt * f.mVeloc + half_dt_dt * f.mAccel;
        
        // zero stiffness for next time step
        f.mStiff.setZero();
    }
}

// update fields on a fluid point
void NewmarkTimeScheme::
update(const std::vector<std::shared_ptr<FluidPoint>> &points) const {
    static const numerical::Real dt = mDt;
    static const numerical::Real half_dt = .5 * mDt;
    static const numerical::Real half_dt_dt = .5 * mDt * mDt;
    for (const std::shared_ptr<FluidPoint> &fp: points) {
        // update dt
        FluidPoint::Fields &f = fp->getFields();
        f.mVeloc += half_dt * (f.mAccel + f.mStiff);
        f.mAccel = f.mStiff;
        f.mDispl += dt * f.mVeloc + half_dt_dt * f.mAccel;
        
        // zero stiffness for next time step
        f.mStiff.setZero();
    }
}
