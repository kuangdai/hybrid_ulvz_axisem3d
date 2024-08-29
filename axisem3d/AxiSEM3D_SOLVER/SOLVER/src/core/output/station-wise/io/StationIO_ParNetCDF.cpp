//
//  StationIO_ParNetCDF.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/24/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  parallel NetCDF IO for station output

#include "StationIO_ParNetCDF.hpp"
#include "io.hpp"
#include "mpi.hpp"
#include "vector_tools.hpp"

// initialize
void StationIO_ParNetCDF::initialize(const std::string &groupName,
                                     int numRecordSteps,
                                     const std::vector<std::string> &channels,
                                     const std::vector<std::string> &stKeys) {
    // finalize
    finalize();
    
    // base
    StationIO::initialize(groupName, numRecordSteps, channels, stKeys);
    if (mNumStationsGlobal == 0) {
        return;
    }
    
    // filename
    const std::string &group_dir = (io::gOutputDirectory + "/stations/"
                                    + groupName);
    const std::string &fname = group_dir + "/axisem3d_synthetics.nc";
    
    // gather stations
    std::vector<std::vector<std::string>> stKeys_rank;
    mpi::gather(stKeys, stKeys_rank, mRankWithMaxNumStations);
    
    /////////////////////// only on max rank ///////////////////////
    if (mpi::rank() == mRankWithMaxNumStations) {
        // flatten keys
        std::vector<std::string> stKeys_all;
        for (int irank = 0; irank < mpi::nproc(); irank++) {
            stKeys_all.insert(stKeys_all.end(),
                              stKeys_rank[irank].begin(),
                              stKeys_rank[irank].end());
        }
        
        // create and open
        mNcFile = std::make_unique<NetCDF_Writer>();
        mNcFile->open(fname, true);
        
        ///////////////////// define variables /////////////////////
        mNcFile->defModeOn();
        
        // time
        mVarID_Time = mNcFile->defineVariable("time_points", {
            numRecordSteps}, numerical::dErr);
        
        // data
        mVarID_Data = mNcFile->defineVariable("data", {
            numRecordSteps, (int)channels.size(), (int)stKeys_all.size()
        }, (numerical::Real)numerical::dErr);
        
        // channels
        int maxChLength = vector_tools::maxLength(channels);
        mNcFile->defineVariable("channel_order", {
            (int)channels.size(), maxChLength}, (char)0);
        
        // stations
        int maxStLength = vector_tools::maxLength(stKeys_all);
        mNcFile->defineVariable("station_order", {
            (int)stKeys_all.size(), maxStLength}, (char)0);
        
        // end defining variables
        mNcFile->defModeOff();
        
        ///////////////////// write info /////////////////////
        // write channels
        for (int ich = 0; ich < channels.size(); ich++) {
            mNcFile->writeVariable("channel_order", channels[ich],
                                   {ich, 0}, {1, (int)channels[ich].size()});
        }
        
        // write station keys
        for (int ist = 0; ist < stKeys_all.size(); ist++) {
            mNcFile->writeVariable("station_order", stKeys_all[ist],
                                   {ist, 0}, {1, (int)stKeys_all[ist].size()});
        }
        
        // close serial
        mNcFile->close();
    }
    
    // open parallel
    mpi::barrier();
    mNcFile->openParallel(fname);
    
    // reset time line
    mFileLineTime = 0;
}

// finalize
void StationIO_ParNetCDF::finalize() {
    // close files
    if (mNcFile) {
        mNcFile->close();
        mNcFile = nullptr;
    }
    mFileLineTime = 0;
}

// dump to file
void StationIO_ParNetCDF::dumpToFile(const eigen::DColX &bufferTime,
                                     const eigen::RTensor3 &bufferFields,
                                     int bufferLine) {
    // no station
    int nst = (int)bufferFields.dimensions()[2];
    if (nst == 0) {
        return;
    }
    
    // no line
    if (bufferLine == 0) {
        return;
    }
    
    // write time
    mNcFile->writeVariable(mVarID_Time, "time_points", bufferTime,
                           {mFileLineTime}, {bufferLine});
    
    // write data
    int nch = (int)bufferFields.dimensions()[1];
    mNcFile->writeVariable(mVarID_Data, "data", bufferFields,
                           {mFileLineTime, 0, mGlobalIndexFirstStation},
                           {bufferLine, nch, nst});
    
    // flush?
    // mNcFile->flush();
    
    // update record postion in file
    mFileLineTime += bufferLine;
}
