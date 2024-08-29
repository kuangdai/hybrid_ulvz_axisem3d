//
//  TimeScheme.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/30/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  time scheme

#ifndef TimeScheme_hpp
#define TimeScheme_hpp

// domain
class SolidPoint;
class FluidPoint;
class Domain;
#include <memory>

// stream
#include "sbstream.hpp"

// verbose
class SimpleTimer;
#include <map>

class TimeScheme {
public:
    // parameter constructor
    TimeScheme(int verboseInterval, int stabilityInterval):
    mVerboseInterval(verboseInterval), mStabilityInterval(stabilityInterval) {
        // nothing
    }
    
    // stream constructor
    TimeScheme(sbs::ifstream &ifs):
    mVerboseInterval(std::min(ifs.get<int>(), 50)), mStabilityInterval(std::min(ifs.get<int>(), 1000)) {
	// nothing
    }
    
    // destructor
    virtual ~TimeScheme() = default;
    
    
    //////////////////// time and domain ////////////////////
    // set time
    void setTime(double t0, double dt, int numTimeSteps) {
        mT0 = t0;
        mDt = dt;
        mNumTimeSteps = numTimeSteps;
    }
    
    // set domain
    void setDomain(const std::shared_ptr<const Domain> &domain) {
        mDomain = domain;
    }
    
    
    //////////////////// point ////////////////////
    // create fields on a solid point
    virtual void createFields(SolidPoint &point) const = 0;
    
    // create fields on a fluid point
    virtual void createFields(FluidPoint &point) const = 0;
    
    
    //////////////////// timeloop ////////////////////
    // solve
    virtual void solve() const = 0;
    
protected:
    // verbose begin
    void verboseBegin(const std::string &ts_name) const;
    // verbose iteration
    void verboseIter(double elapsed, int tstep, double t) const;
    // verbose end
    void verboseEnd(const std::string &ts_name,
                    double elapsed, double t) const;
    
    
    //////////////////// data ////////////////////
    // intervals
    const int mVerboseInterval;
    const int mStabilityInterval;
    
    // time
    double mT0 = 0.;
    double mDt = 0.;
    int mNumTimeSteps = 0;
    
    // domain
    std::shared_ptr<const Domain> mDomain = nullptr;
    
    
    ////////////////////////////////////////
    //////////////// static ////////////////
    ////////////////////////////////////////
    
    // report loop timers
    static
    void reportLoopTimers(const std::vector<std::string> &timerNames,
                          const std::map<std::string, SimpleTimer> &timers);
    
public:
    // stream-based creator
    static std::unique_ptr<TimeScheme> createFromStream(sbs::ifstream &ifs);
};

#endif /* TimeScheme_hpp */
