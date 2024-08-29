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
   

    if (z < 0) return 1;
    return 810;


 
    // output: Fourier Expansion order nu
    int nu = 2;
    
    ////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////
    // TODO: Compute nu from s, z, r, theta, depth and mParameters
    // NOTE: Only edit within this box!
    //       If this is the first time you are looking at this part
    //       of code, the following example implements a Nu field
    //       that only depends on radius (depth), with a very large Nu
    //       assigned near the core-mantle boundary
    
    // double rcmb = 3480e3;
    // double r_low_broad = rcmb - 500e3;
    // double r_low_narrow = rcmb - 100e3;
    // double r_upp_narrow = rcmb + 200e3;
    // double r_upp_broad = rcmb + 1000e3;
    // int nu_narrow = 1000;
    // int nu_broad = 100; 
    double r_earth = 6371e3;
    double r_high = r_earth - 400e3;
    // I used 400 for the benchmark of noise-correlations with specfem, long period
    double r_low = r_earth - 600e3;
    // I used 600 km for the long-period specfem noise-corr benchmark
    int nu_high = 300;
    // I used 300 for the longperiod specfem noise-corr bm with s40rts
    int nu_low = 0;

    // first, set low nu at ep. distance larger 90 degree (laura)
    if (theta > pi / 4) return 1;
    // for the specfem noise-correlation benchmark, 45 degrees was used. (pi / 4)

    if (r > r_high) {
        nu = nu_high;
        // std::cout << "above r_high" << std::endl;
    } else if (r < r_low) {
        nu = nu_low;
      //  std::cout << "below r_low" << std::endl;
    } else {
       // std::cout << "in between" << std::endl;
        nu = (nu_high - nu_low) / (r_high - r_low) * (r - r_low) + nu_low;
    }

    // Kuangdai's ULVZ nu:
    // if (r <= r_low_broad) {
    //    nu = nu_broad;
    // } else if (r <= r_low_narrow) {
    //    nu = (nu_narrow - nu_broad) / (r_low_narrow - r_low_broad) * (r - r_low_broad) + nu_broad;
    // } else if (r <= r_upp_narrow) {
    //    nu = nu_narrow;
    // } else if (r <= r_upp_broad) {
    //    nu = (nu_narrow - nu_broad) / (r_upp_narrow - r_upp_broad) * (r - r_upp_broad) + nu_broad;
    // } else {
    //    nu = nu_broad;
    // }
    
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

