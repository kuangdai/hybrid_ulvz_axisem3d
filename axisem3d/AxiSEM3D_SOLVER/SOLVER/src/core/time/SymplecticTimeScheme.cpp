//
//  SymplecticTimeScheme.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/30/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  symplectic time scheme

#include "SymplecticTimeScheme.hpp"
#include "SolidPoint.hpp"
#include "FluidPoint.hpp"
#include "timer.hpp"
#include "Domain.hpp"

//////////////////// point ////////////////////
// create fields on a solid point
void SymplecticTimeScheme::createFields(SolidPoint &point) const {
    int nu_1 = point.getNu() + 1;
    SolidPoint::Fields &f = point.getFields();
    f.mStiff = eigen::CMatX3::Zero(nu_1, 3);
    f.mDispl = eigen::CMatX3::Zero(nu_1, 3);
    f.mVeloc = eigen::CMatX3::Zero(nu_1, 3);
}

// create fields on a fluid point
void SymplecticTimeScheme::createFields(FluidPoint &point) const {
    int nu_1 = point.getNu() + 1;
    FluidPoint::Fields &f = point.getFields();
    f.mStiff = eigen::CColX::Zero(nu_1, 1);
    f.mDispl = eigen::CColX::Zero(nu_1, 1);
    f.mVeloc = eigen::CColX::Zero(nu_1, 1);
}


//////////////////// timeloop ////////////////////
// solve
void SymplecticTimeScheme::solve() const {
    // loop timer
    std::vector<std::string> timerNames = {"TOTAL", "SOURCE",
        "STIFFNESS", "MASS_TERM", "TIME_MARCH", "SF_COUPLING",
        "BOUNDARIES", "MPI_COMM", "WAVE_OUTPUT"};
    std::map<std::string, SimpleTimer> timers;
    for (const std::string &name: timerNames) {
        timers.insert({name, SimpleTimer()});
    }
    
    // verbose
    verboseBegin("SYMPLECTIC");
    
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
        // wavefield output: record and dump
        // Need to do this here to have the same time-axis with Newmark
        timers.at("WAVE_OUTPUT").resume();
        mDomain->recordStations(tstep, t);
        timers.at("WAVE_OUTPUT").pause();
        
        // apply source
        timers.at("SOURCE").resume();
        mDomain->applySources(tstep, t);
        timers.at("SOURCE").pause();
        
        // launch
        timers.at("TIME_MARCH").resume();
        launch(mDomain->getSolidPoints());
        launch(mDomain->getFluidPoints());
        timers.at("TIME_MARCH").pause();
        
        // sub iteration
        for (int isub = 1; isub <= 4; isub++) {
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
            
            // stability
            if (tstep % mStabilityInterval == 0 && isub == 1) {
                mDomain->checkStability(tstep + 1, t, mDt);
            }
            
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
            update(mDomain->getSolidPoints(), isub);
            update(mDomain->getFluidPoints(), isub);
            timers.at("TIME_MARCH").pause();
        } // end of sub iteration
        
        // verbose
        verboseIter(timers.at("TOTAL").elapsedTotal(), tstep + 1, t);
        
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
    verboseEnd("SYMPLECTIC", timers.at("TOTAL").elapsedTotal(), t);
    
    // report timers
    reportLoopTimers(timerNames, timers);
}

// launch on solid points
void SymplecticTimeScheme::
launch(const std::vector<std::shared_ptr<SolidPoint>> &points) const {
    static const numerical::Real kappa_dt = sKappa[0] * mDt;
    for (const std::shared_ptr<SolidPoint> &sp: points) {
        SolidPoint::Fields &f = sp->getFields();
        f.mDispl += kappa_dt * f.mVeloc;
    }
}

// launch on fluid points
void SymplecticTimeScheme::
launch(const std::vector<std::shared_ptr<FluidPoint>> &points) const {
    static const numerical::Real kappa_dt = sKappa[0] * mDt;
    for (const std::shared_ptr<FluidPoint> &fp: points) {
        FluidPoint::Fields &f = fp->getFields();
        f.mDispl += kappa_dt * f.mVeloc;
    }
}

// update fields on solid points
void SymplecticTimeScheme::
update(const std::vector<std::shared_ptr<SolidPoint>> &points,
       int iSubIter) const {
    const numerical::Real pi_dt = sPi[iSubIter] * mDt;
    const numerical::Real kappa_dt = sKappa[iSubIter] * mDt;
    for (const std::shared_ptr<SolidPoint> &sp: points) {
        SolidPoint::Fields &f = sp->getFields();
        // update dt
        f.mVeloc += pi_dt * f.mStiff;
        f.mDispl += kappa_dt * f.mVeloc;
        
        // zero stiffness for next time step
        f.mStiff.setZero();
    }
}

// update fields on fluid points
void SymplecticTimeScheme::
update(const std::vector<std::shared_ptr<FluidPoint>> &points,
       int iSubIter) const {
    const numerical::Real pi_dt = sPi[iSubIter] * mDt;
    const numerical::Real kappa_dt = sKappa[iSubIter] * mDt;
    for (const std::shared_ptr<FluidPoint> &fp: points) {
        FluidPoint::Fields &f = fp->getFields();
        // update dt
        f.mVeloc += pi_dt * f.mStiff;
        f.mDispl += kappa_dt * f.mVeloc;
        
        // zero stiffness for next time step
        f.mStiff.setZero();
    }
}
