//
//  io.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/8/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  IO interfaces

#include "io.hpp"
#include "bstring.hpp"

// io::cout
#include <iostream>

// mkdir
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
};
#include "mpi.hpp"

namespace io {
    //////////////////// input/output dirs ////////////////////
    std::string gInputDirectory = "";
    std::string gOutputDirectory = "";
    
    // check dir existence
    bool dirExists(const std::string &path) {
        struct stat info;
        if (stat(path.c_str(), &info) != 0) {
            return false;
        } else if (info.st_mode & S_IFDIR) {
            return true;
        } else {
            return false;
        }
    }
    
    // mkdir
    void mkdir(const std::string &path) {
        if (!dirExists(path)) {
            ::mkdir(path.c_str(), ACCESSPERMS);
        }
    }
    
    // verify input/output dirs under executable dir
    void verifyDirectories(int argc, char *argv[], std::string &warning) {
        // find path of executable
        std::string argv0(argv[0]);
        std::size_t found = argv0.find_last_of("/\\");
        if (found == std::string::npos) {
            argv0 = "./" + argv0;
            found = argv0.find_last_of("/\\");
        }
        std::string execDirectory = argv0.substr(0, found);
        
        // so far, this problem happens with valgrind only
        if (execDirectory.size() == 0) {
            execDirectory = ".";
        }
        
        // input and output
        gInputDirectory = execDirectory + "/input";
        gOutputDirectory = execDirectory + "/output";
        warning = "";
        if (mpi::root()) {
            if (!dirExists(gInputDirectory)) {
                throw std::runtime_error("io::verifyDirectories || "
                                         "Missing input directory: ||"
                                         + gInputDirectory);
            }
            if (dirExists(gOutputDirectory)) {
                // backup the old? maybe unnecessary, use C rename() if needed
                warning = bstring::warning("io::verifyDirectories || "
                                           "Output directory exists: ||" +
                                           gOutputDirectory + "|| "
                                           "Old results will be overwritten.");
            }
            mkdir(gOutputDirectory);
            mkdir(gOutputDirectory + "/stations");
            mkdir(gOutputDirectory + "/develop");
        }
    }
    
    
    //////////////////// runtime verbose ////////////////////
    // verbose control
    VerboseLevel gVerbose = VerboseLevel::detailed;
    bool gVerboseWarnings = true;
    
    // print welcome messages
    std::string welcome() {
        std::string welc = "\n"
        "{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~}"
        "[                                                            ]"
        "[    A            |i .|'''.|'||''''E'||    ||'  ____'||''|.  ]"
        "[   |||   ... ...... ||..  ' ||  .   |||  |||   ` // ||   || ]"
        "[  |  ||   '|..'  ||  ''|||. ||''|   |'|..'||    //  ||    ||]"
        "[ .''''|.   .x.   ||      '||||      | 'M' ||    \\  ||    ||]"
        "[.|.  .||..|  ||..||.|'...|S.||....|.|. | .||.    3' D|...|' ]"
        "[.............................................   //          ]"
        "[                                               /'     v x.y ]"
        "[                                                            ]"
        "[Copyright (c) 2019 Kuangdai Leng & friends, MIT License     ]"
        "[News, suggestions and issues: www.axisem3d.ox.ac.uk         ]"
        "{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~}\n\n";
        std::string space = bstring::filled((bstring::fmt::gBoxWidth - 60) / 2);
        std::string tilde = bstring::filled((bstring::fmt::gBoxWidth - 60) / 2,
                                            '~');
        welc = bstring::replace(welc, "{", tilde);
        welc = bstring::replace(welc, "}", tilde + "\n");
        welc = bstring::replace(welc, "[", space);
        welc = bstring::replace(welc, "]", space + "\n");
        welc = bstring::replace(welc, "x.y", _VERSION);
        return bstring::replace(welc, "\\", "\\\\");
    }
    
    // verbose
    std::string verbose() {
        std::stringstream ss;
        ss << bstring::boxTitle("IO");
        ss << bstring::boxSubTitle(0, "Directories");
        ss << bstring::boxEquals(2, 8, "input", gInputDirectory);
        ss << bstring::boxEquals(2, 8, "output", gOutputDirectory);
        ss << bstring::boxEquals(2, 8, "source", gProjectDirectory);
        ss << bstring::boxSubTitle(0, "Verbose");
        if (gVerbose == VerboseLevel::essential) {
            ss << bstring::boxEquals(2, 8, "level", "essential");
        } else if (gVerbose == VerboseLevel::detailed) {
            ss << bstring::boxEquals(2, 8, "level", "detailed");
        } else {
            ss << bstring::boxEquals(2, 8, "level", "none");
        }
        ss << bstring::boxEquals(2, 8, "warnings", gVerboseWarnings);
        ss << bstring::boxBaseline() << "\n\n";
        return ss.str();
    }
    
    ////////////////////////////// cout on root //////////////////////////////
    mpi_root_cout::mpi_root_cout()
    :mMyWorldRank(-1), mCoutWorldRank(0), mCoutStream(&(std::cout)) {
        // nothing
    }
    mpi_root_cout cout;
}
