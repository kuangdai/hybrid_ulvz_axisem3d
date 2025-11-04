//
//  TimeScheme.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/2/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

#include "TimeScheme.hpp"
// verbose
#include "io.hpp"
#include "bstring.hpp"
// timer
#include "timer.hpp"
#include "eigen_generic.hpp"
#include "mpi.hpp"
#include <iomanip>
// stream
#include "NewmarkTimeScheme.hpp"
#include "SymplecticTimeScheme.hpp"

using namespace bstring;

// verbose begin
void TimeScheme::verboseBegin(const std::string &ts_name) const {
    if (io::gVerbose != io::VerboseLevel::none) {
        std::stringstream ss;
        ss << io::endl << boxBaseline(fmt::gBoxWidth, 'T');
        ss << boxTitle(ts_name + " TIME LOOP STARTS", fmt::gBoxWidth, 'T');
        ss << boxBaseline(fmt::gBoxWidth, 'T') << io::endl << io::endl;
        io::cout << ss.str();
    }
}

// verbose iteration
void TimeScheme::verboseIter(double elapsed, int tstep, double t) const {
    if (tstep % mVerboseInterval == 0 &&
        io::gVerbose != io::VerboseLevel::none) {
        elapsed /= 3600.;
        double total = elapsed / tstep * mNumTimeSteps;
        double left = total - elapsed;
        int percent = (int)(100. * tstep / mNumTimeSteps);
        std::stringstream ss;
        std::stringstream sst;
        sst << tstep << " / " << mNumTimeSteps;
        sst << " (" << percent << "%)";
        ss << boxEquals(2, 27, "wave propagation time / sec", t);
        ss << boxEquals(2, 27, "time step / num total steps", sst.str());
        ss << boxEquals(2, 27, "elapsed wall-clock time / h", elapsed);
        ss << boxEquals(2, 27, "est. remaining w-c time / h", left);
        ss << boxEquals(2, 27, "est. aggregate w-c time / h", total);
        ss << io::endl;
        io::cout << ss.str();
    }
}

// verbose end
void TimeScheme::verboseEnd(const std::string &ts_name,
                            double elapsed, double t) const {
    if (io::gVerbose != io::VerboseLevel::none) {
        std::stringstream ss;
        ss << io::endl << boxBaseline(fmt::gBoxWidth, 'T');
        ss << boxTitle(ts_name + " TIME LOOP FINISHES", fmt::gBoxWidth, 'T');
        ss << boxBaseline(fmt::gBoxWidth, 'T') << io::endl;
        ss << boxEquals(0, 27, "wave propagation time / sec", t);
        ss << boxEquals(0, 27, "total # time steps finished", mNumTimeSteps);
        ss << boxEquals(0, 27, "elapsed wall-clock time / h", elapsed / 3600.);
        ss << io::endl;
        io::cout << ss.str();
    }
}


////////////////////////////////////////
//////////////// static ////////////////
////////////////////////////////////////

// report loop timers
void TimeScheme::
reportLoopTimers(const std::vector<std::string> &timerNames,
                 const std::map<std::string, SimpleTimer> &timers) {
    // collect local timers
    eigen::DMatXX elapsed(mpi::nproc(), timerNames.size());
    elapsed.setZero();
    for (int itimer = 0; itimer < timerNames.size(); itimer++) {
        elapsed(mpi::rank(), itimer) =
        timers.at(timerNames[itimer]).elapsedTotal();
    }
    
    // collect timers over ranks
    mpi::sumEigen(elapsed);
    
    // write to file
    if (mpi::root()) {
        std::string fname = io::gOutputDirectory + "/develop/loop_timers.log";
        std::ofstream ofs(fname);
        if (!ofs) {
            throw std::runtime_error("TimeScheme::reportLoopTimers || "
                                     "Error creating log file: || " + fname);
        }
        // header
        int width = 15;
        ofs << std::setw(width) << "MPI_RANK";
        for (int itimer = 0; itimer < timerNames.size(); itimer++) {
            ofs << std::setw(width) << timerNames[itimer];
        }
        ofs << io::endl;
        // data: eigen does not support width
        for (int row = 0; row < elapsed.rows(); row++) {
            ofs << std::setw(width) << row; // mpi rank
            for (int col = 0; col < elapsed.cols(); col++) {
                ofs << std::setw(width) << elapsed(row, col);
            }
            ofs << io::endl;
        }
        ofs << io::endl;
    }
}

// stream-based creator
std::unique_ptr<TimeScheme> TimeScheme::createFromStream(sbs::ifstream &ifs) {
    const std::string &clsName = ifs.get<std::string>();
    if (clsName == "NewmarkTimeScheme") {
        return std::make_unique<NewmarkTimeScheme>(ifs);
    } else if (clsName == "SymplecticTimeScheme") {
        return std::make_unique<SymplecticTimeScheme>(ifs);
    } else {
        throw std::runtime_error("TimeScheme::createFromStream || "
                                 "Unknown derived class of TimeScheme: "
                                 + clsName);
    }
}
