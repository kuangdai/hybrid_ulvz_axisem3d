//
//  main.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/6/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  AxiSEM3D main

#include "mpi.hpp"
#include "io.hpp"
#include "inparam.hpp"
#include "timer.hpp"
#include "bstring.hpp"
#include <iostream>

#include "ExodusMesh.hpp"
#include "geodesy.hpp"


extern "C" void set_ftz();
extern void main_stream_1_0(int argc, char *argv[]);

int main(int argc, char *argv[]) {
    // denormal float handling
    set_ftz();
    
    // initialise MPI
    mpi::initialize(argc, argv);
    
    // all other codes go into try
    try {
        // initialise IO
        std::string warningIO;
        io::verifyDirectories(argc, argv, warningIO);
        
        // inparam setup
        inparam::setup();
        
        // welcome
        if (io::gVerbose != io::VerboseLevel::none) {
            io::cout << io::welcome();
        }
        
        // verbose mpi, io, inparam
        if (io::gVerbose == io::VerboseLevel::detailed) {
            io::cout << mpi::verbose();
            io::cout << io::verbose();
            io::cout << inparam::verbose();
        }
        
        // io warning
        if (io::gVerboseWarnings) {
            io::cout << warningIO;
        }
        
        // preloop timer
        timer::setup(inparam::gInparamAdvanced.gets<bool>
                     ("develop:diagnose_preloop"));
        if (io::gVerbose == io::VerboseLevel::detailed) {
            io::cout << timer::verbose();
        }
        
        // exodus mesh and geodesy
        ExodusMesh exodusMesh(inparam::gInparamModel);
        geodesy::setup(inparam::gInparamModel, exodusMesh);
        if (io::gVerbose != io::VerboseLevel::none) {
            io::cout << exodusMesh.verbose();
            io::cout << geodesy::verbose(exodusMesh.getDiscontinuities());
        }
        
        // run from stream created by 1.0
        main_stream_1_0(argc, argv);
    } catch (const std::runtime_error &e) {
        std::cout << bstring::exception(e);
        mpi::abort(1);
    }
    
    // finalize mpi
    mpi::finalize();
    return 0;
}
