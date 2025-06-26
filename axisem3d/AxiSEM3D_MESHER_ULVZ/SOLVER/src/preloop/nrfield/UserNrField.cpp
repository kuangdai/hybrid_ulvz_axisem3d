// UserNrField.cpp
// created by Kuangdai on 11-Jun-2017 
// user-defined nr integer field

#include <iostream>
#include "UserNrField.h"
#include <sstream>
#include "Geodesy.h"

UserNrField::UserNrField(bool useLucky, const std::vector<double> &params): 
NrField(useLucky), mParameters(params) {
    // nothing
}

int UserNrField::getNrAtPoint(const RDCol2 &coords) const {
    // input: coordinates
    // s, z, r, theta, depth
    double s = coords(0);
    double z = coords(1);
    double r, theta;
    Geodesy::rtheta(coords, r, theta);
    double depth = Geodesy::getROuter() - r;
   
    // output: Fourier Expansion order nu
    int nu = mParameters[0];
    
    ////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////
    // TODO: Compute nu from s, z, r, theta, depth and mParameters
    // NOTE: Only edit within this box!
    //       If this is the first time you are looking at this part
    //       of code, the following example implements a Nu field
    //       that only depends on radius (depth), with a very large Nu
    //       assigned near the core-mantle boundary
    
    int nu_big = mParameters[0];
    int nu_small = 5;
    if (mParameters.size() >= 2) {
        nu_small = mParameters[1];
    }
    double rcmb = 3480e3;
    double r_low0 = rcmb - 100e3;
    if (mParameters.size() >= 3) {
        r_low0 = mParameters[2];
    }
    double r_low1 = rcmb - 300e3;
    if (mParameters.size() >= 4) {
        r_low1 = mParameters[3];
    }
    double r_upp0 = rcmb + 200e3;
    if (mParameters.size() >= 5) {
        r_upp0 = mParameters[4];
    }
    double r_upp1 = rcmb + 500e3;
    if (mParameters.size() >= 6) {
        r_upp1 = mParameters[5];
    }
    double theta0 = 0.25;
    if (mParameters.size() >= 7) {
        theta0 = mParameters[6];
    }
    double theta1 = 0.5;
    if (mParameters.size() >= 8) {
        theta1 = mParameters[7];
    }

    if (r >= r_low0 && r <= r_upp0 && theta <= theta0) {
        nu = nu_big;
    } else if (r <= r_low1 || r >= r_upp1 || theta >= theta1) {
        nu = nu_small;
    } else {
        // transition area
        // interpolate based on radius and theta
        double fr = 0.0;
        if (r < r_low0) {
            fr = (r - r_low1) / (r_low0 - r_low1); // from 0 to 1
        } else if (r > r_upp0) {
            fr = (r_upp1 - r) / (r_upp1 - r_upp0); // from 1 to 0
        } else {
            fr = 1.0;
        }

        double ft = 1.0;
        if (theta > theta0 && theta < theta1) {
            ft = (theta1 - theta) / (theta1 - theta0);
        }

        double weight = fr * ft;
        nu = (int)(nu_small + weight * (nu_big - nu_small));
    }
    // NOTE: Only edit within this box!
    ////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    // get nr from nu
    int nr = 2 * nu + 1;
    if (nr <= 0) {
        throw std::runtime_error("UserNrField::getNrAtPoint || Non-positive Nr.");
    }
    return nr;
}

std::string UserNrField::verbose() const {
    std::stringstream ss;
    ss << "\n================= Fourier Expansion Order ==================" << std::endl;
    ss << "  Type                     =   User-defined" << std::endl;
    if (mParameters.size() > 0) {
        ss << "  Parameter List           =   ";
        for (int ipar = 0; ipar < mParameters.size(); ipar++) {
            ss << mParameters[ipar] << " ";
        }
        ss << std::endl;
    }
    ss << "  Use FFTW Lucky Numbers   =   " << (mUseLuckyNumber ? "YES" : "NO") << std::endl;
    ss << "================= Fourier Expansion Order ==================\n" << std::endl;
    return ss.str();
}

