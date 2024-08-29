//
//  small_tests.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/13/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  add any lightweight test functions here

#include "eigen_tools.hpp"
#include "eigen_generic.hpp"
#include "eigen_point.hpp"
#include "eigen_element.hpp"
#include "timer.hpp"
#include "fftw.hpp"
#include <iostream>

using numerical::Real;
using numerical::dPi;

void test_fft() {
    // fftw
    int nrOdd = 5;
    int nrEven = nrOdd - 1;
    int ncOdd = nrOdd / 2 + 1;
    int ncEven = nrEven / 2 + 1;
    fftw::gFFT_1.addNR(nrOdd);
    fftw::gFFT_1.addNR(nrEven);
    fftw::gFFT_1.addNR(nrOdd - 2);
    fftw::gFFT_1.createPlans(1.);
    
    // test fft odd
    eigen::CColX cOdd = eigen::CColX::Random(ncOdd);
    cOdd.row(0).imag().setZero();
    eigen::CColX cOdd_back = eigen::CColX::Random(ncOdd);
    eigen::RColX rOdd = eigen::RColX::Zero(nrOdd);
    fftw::gFFT_1.computeC2R(cOdd, rOdd, nrOdd);
    fftw::gFFT_1.computeR2C(rOdd, cOdd_back, nrOdd);
    std::cout << cOdd.transpose() << std::endl;
    std::cout << cOdd_back.transpose() << std::endl << std::endl;
    
    // test fft even
    eigen::CColX cEven = eigen::CColX::Random(ncEven);
    cEven.row(0).imag().setZero();
    cEven.row(ncEven - 1).imag().setZero();
    eigen::CColX cEven_back = eigen::CColX::Random(ncEven);
    eigen::RColX rEven = eigen::RColX::Zero(nrEven);
    fftw::gFFT_1.computeC2R(cEven, rEven, nrEven);
    fftw::gFFT_1.computeR2C(rEven, cEven_back, nrEven);
    std::cout << cEven.transpose() << std::endl;
    std::cout << cEven_back.transpose() << std::endl << std::endl;
    
    // test eigen_tools odd
    typedef Eigen::Matrix<Real, 1, 1> RRow1;
    for (int ir = 0; ir < nrOdd; ir++) {
        RRow1 valueAtPhi;
        Real dphi = 2. * dPi / nrOdd;
        eigen_tools::computeFourierAtPhi<Real, 1>(cOdd_back, ncOdd, ir * dphi,
                                                  valueAtPhi);
        std::cout << valueAtPhi(0) << " ";
    }
    std::cout << std::endl;
    std::cout << rOdd.transpose() << std::endl << std::endl;
    
    // test eigen_tools even
    typedef Eigen::Matrix<Real, 1, 1> RRow1;
    for (int ir = 0; ir < nrEven; ir++) {
        RRow1 valueAtPhi;
        Real dphi = 2. * dPi / nrEven;
        eigen_tools::computeFourierAtPhi<Real, 1>(cEven_back, ncEven, ir * dphi,
                                                  valueAtPhi);
        std::cout << valueAtPhi(0) << " ";
    }
    std::cout << std::endl;
    std::cout << rEven.transpose() << std::endl << std::endl;
    
    // test odd -> odd
    int ncOddS = ncOdd - 1;
    int nrOddS = nrOdd - 2;
    eigen::CColX cOddS = cOdd.topRows(ncOddS);
    eigen::RColX rOddS = eigen::RColX::Zero(nrOddS);
    for (int ir = 0; ir < nrOddS; ir++) {
        RRow1 valueAtPhi;
        Real dphi = 2. * dPi / nrOddS;
        eigen_tools::computeFourierAtPhi<Real, 1>(cOddS, ncOddS, ir * dphi,
                                                  valueAtPhi);
        std::cout << valueAtPhi(0) << " ";
    }
    std::cout << std::endl;
    fftw::gFFT_1.computeC2R(cOddS, rOddS, nrOddS);
    std::cout << rOddS.transpose() << std::endl << std::endl;
    
    // test odd -> even
    int ncEvenS = ncOdd;
    int nrEvenS = nrOdd - 1;
    eigen::CColX cEvenS = cOdd.topRows(ncEvenS);
    cEvenS.bottomRows(1).imag().setZero();
    eigen::RColX rEvenS = eigen::RColX::Zero(nrEvenS);
    for (int ir = 0; ir < nrEvenS; ir++) {
        RRow1 valueAtPhi;
        Real dphi = 2. * dPi / nrEvenS;
        eigen_tools::computeFourierAtPhi<Real, 1>(cOdd, ncOdd, ir * dphi,
                                                  valueAtPhi);
        std::cout << valueAtPhi(0) << " ";
    }
    std::cout << std::endl;
    fftw::gFFT_1.computeC2R(cEvenS, rEvenS, nrEvenS);
    std::cout << rEvenS.transpose() << std::endl << std::endl;
}

void performance_fft() {
    SimpleTimer timer;
    
    // fftw
    int nr = 9;
    int nc = nr / 2 + 1;
    fftw::gFFT_N3.addNR(nr);
    fftw::gFFT_N3.createPlans(1.);
    
    // test fft odd
    eigen::CMatXN3 cdata = eigen::CMatXN3::Random(nc, 3 * spectral::nPEM);
    eigen::RMatXN3 rdata = eigen::RMatXN3::Zero(nr, 3 * spectral::nPEM);
    
    int N = 1000000;
    timer.start();
    for (int i = 0; i < N; i++) {
        fftw::gFFT_N3.computeC2R(cdata, rdata, nr);
        fftw::gFFT_N3.computeR2C(rdata, cdata, nr);
    }
    timer.pause();
    std::cout << timer.elapsedTotal() << std::endl;
}

void performance_eigen() {
    SimpleTimer timer;
    
    // R = F * P, with P in diffferent structures
    eigen::CMatPP_RM F, R;
    eigen::RMatPP_RM P;
    Eigen::Matrix<numerical::Real, 1, spectral::nPEM> P1;
    F.setRandom();
    P.setRandom();
    P1.setRandom();
    
    int N = 1000000;
    
    // P in 5 * 5
    timer.start();
    for (int i = 0; i < N; i++) {
        R = F.cwiseProduct(P);
    }
    timer.pause();
    std::cout << timer.elapsedTotal() << std::endl;
    
    // P in 1 * 25
    timer.start();
    for (int i = 0; i < N; i++) {
        R = F.cwiseProduct(Eigen::Map<eigen::RMatPP_RM>(P1.data()));
    }
    timer.pause();
    std::cout << timer.elapsedTotal() << std::endl;
}
