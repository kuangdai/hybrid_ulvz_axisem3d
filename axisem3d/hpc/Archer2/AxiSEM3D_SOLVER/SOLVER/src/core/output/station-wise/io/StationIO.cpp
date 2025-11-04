//
//  StationIO.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/24/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  IO for station output

#include "StationIO.hpp"
#include "io.hpp"
#include "mpi.hpp"
#include <numeric>
#include <fstream>

void StationIO::initialize(const std::string &groupName,
                           int numRecordSteps,
                           const std::vector<std::string> &channels,
                           const std::vector<std::string> &stKeys) {
    // mpi globals
    int nst = (int)stKeys.size();
    std::vector<int> nstg;
    mpi::gather(nst, nstg, -1);
    mNumStationsGlobal = std::accumulate(nstg.begin(), nstg.end(), 0);
    if (mNumStationsGlobal == 0) {
        // no station at all
        mRankWithMaxNumStations = -1;
        mGlobalIndexFirstStation = 0;
    } else {
        mRankWithMaxNumStations =
        (int)(std::max_element(nstg.begin(), nstg.end()) - nstg.begin());
        mGlobalIndexFirstStation =
        std::accumulate(nstg.begin(), nstg.begin() + mpi::rank(), 0);
    }
    
    // create group dir
    if (mpi::root()) {
        const std::string &group_dir = (io::gOutputDirectory + "/stations/"
                                        + groupName);
        io::mkdir(group_dir);
    }
    // wait for mkdir
    mpi::barrier();
}

// create rank_station.info
void StationIO::
createRankStation(const std::string &groupName,
                  const std::vector<std::string> &stKeys) const {
    std::vector<std::vector<std::string>> stKeysG;
    mpi::gather(stKeys, stKeysG, 0);
    if (mpi::root()) {
        const std::string &group_dir = (io::gOutputDirectory + "/stations/"
                                        + groupName);
        std::ofstream fout(group_dir + "/rank_station.info");
        if (!fout) {
            throw std::runtime_error("StationIO::initialize || "
                                     "Error creating rank-station index file:"
                                     " || " + group_dir + "/rank_station.info");
        }
        fout << "MPI_RANK STATION_KEY STATION_INDEX_IN_RANK" << "\n";
        for (int rank = 0; rank < mpi::nproc(); rank++) {
            for (int ist = 0; ist < stKeysG[rank].size(); ist++) {
                fout << rank << " " << stKeysG[rank][ist] << " " << ist << "\n";
            }
        }
    }
}
