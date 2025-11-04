//
//  SymplecticTimeScheme.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/30/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  symplectic time scheme

#ifndef SymplecticTimeScheme_hpp
#define SymplecticTimeScheme_hpp

#include "TimeScheme.hpp"
#include "numerical.hpp"

class SymplecticTimeScheme: public TimeScheme {
public:
    // parameter constructor
    SymplecticTimeScheme(int verboseInterval, int stabilityInterval):
    TimeScheme(verboseInterval, stabilityInterval) {
        // nothing
    }
    
    // stream constructor
    SymplecticTimeScheme(sbs::ifstream &ifs):
    TimeScheme(ifs) {
        // nothing
    }
    
    
    //////////////////// point ////////////////////
    // create fields on a solid point
    void createFields(SolidPoint &point) const;
    
    // create fields on a fluid point
    void createFields(FluidPoint &point) const;
    
    
    //////////////////// timeloop ////////////////////
    // solve
    void solve() const;
    
private:
    // launch on solid points
    void launch(const std::vector<std::shared_ptr<SolidPoint>> &points) const;
    
    // launch on fluid points
    void launch(const std::vector<std::shared_ptr<FluidPoint>> &points) const;
    
    // update fields on solid points
    void update(const std::vector<std::shared_ptr<SolidPoint>> &points,
                int iSubIter) const;
    
    // update fields on fluid points
    void update(const std::vector<std::shared_ptr<FluidPoint>> &points,
                int iSubIter) const;
    
    
    //////////////////////////////////////////
    ////////////////// static ////////////////
    //////////////////////////////////////////
    
private:
    inline static const double sAlpha = 0.1786178958448091;
    inline static const double sBeta = -0.2123418310626054;
    inline static const double sGamma = -0.06626458266981849;
    
    inline static const std::vector<numerical::Real> sKappa = {
        (numerical::Real)sAlpha,
        (numerical::Real)sGamma,
        (numerical::Real)(1. - 2. * (sGamma + sAlpha)),
        (numerical::Real)sGamma,
        (numerical::Real)sAlpha};
    
    inline static const std::vector<numerical::Real> sPi = {
        (numerical::Real)0., // unused front padding
        (numerical::Real)(0.5 - sBeta),
        (numerical::Real)sBeta,
        (numerical::Real)sBeta,
        (numerical::Real)(0.5 - sBeta)};
};

#endif /* SymplecticTimeScheme_hpp */
