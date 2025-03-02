//
//  main_stream_1_0.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/12/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  run simulation from bin streams created by AxiSEM3D 1.0

#include "mpi.hpp"
#include "io.hpp"
#include "sbstream.hpp"
#include <fstream>

#include "Domain.hpp"
#include "TimeScheme.hpp"
#include "SourceTimeFunction.hpp"
#include "SolidForce.hpp"
#include "WavefieldInjection.hpp"
#include "StationIO_Ascii.hpp"
#include "StationIO_NetCDF.hpp"
#include "StationIO_ParNetCDF.hpp"
#include "SolidElement.hpp"
#include "StationSolid.hpp"
#include "StationGroup.hpp"

#include "fftw.hpp"
#include "GradientQuadrature.hpp"
#include "se_constants.hpp"
#include <iostream>

#include <cmath>

void main_stream_1_0(int argc, char *argv[]) {
    // input files
    std::ifstream pfs(io::gInputDirectory + "/inparam.stream");
    std::string prompt;
    std::string binDir;
    pfs >> prompt >> binDir;
    binDir = io::gInputDirectory + "/" + binDir;
    sbs::ifstream ifs(binDir + "/rank" + mpi::strRank() + ".bin");
    
    /////////////////// domain ///////////////////
    // time scheme
    auto timeScheme = TimeScheme::createFromStream(ifs);

    // domain
    auto domain = Domain::createFromStream(ifs, *timeScheme);
    if (io::gVerbose != io::VerboseLevel::none) {
        io::cout << domain->verbose();
    }
    
    // set domain
    timeScheme->setDomain(domain);
    
    
    // MESHING STAGE //
    ///////////////////
    // SOLVING STAGE //
    
    
    /////////////////// injection ///////////////////
    std::unique_ptr<WavefieldInjection> wj = nullptr;
    // has injection
    bool HAS_INJECT;
    pfs >> prompt >> HAS_INJECT;
    if (HAS_INJECT) {
        //////////////// init injection ////////////////
        // injeciton par
        std::ifstream wjfs(io::gInputDirectory + "/inparam.stream_wj");
        // create
        wj = std::make_unique<WavefieldInjection>();
        // number of quad
        int nInQuad;
        wjfs >> prompt >> nInQuad;
        // quad tag
        std::vector<int> inQuad(nInQuad, -1);
        wjfs >> prompt;
        for (int iq = 0; iq < nInQuad; iq++) {
            wjfs >> inQuad[iq];
        }
        // boundary location
        double r0, r1, t1, dtol;
        wjfs >> prompt >> r0;
        wjfs >> prompt >> r1;
        wjfs >> prompt >> t1;
        wjfs >> prompt >> dtol;
        wj->initializeElements(inQuad, {r0, r1}, {t1}, dtol, *domain);
       

        // write info
        io::mkdir(io::gOutputDirectory + "/WJ_BOX");
        mpi::barrier();
        std::vector<std::string> quadKeys;
        std::vector<int> quadNrs;
        std::vector<std::vector<int>> quadBoundaryPnts;
        std::vector<std::vector<std::string>> recKeys;
        std::vector<eigen::DMatXX> rec_spz;
        
        // solid
        wj->infoSolid(quadKeys, quadNrs, quadBoundaryPnts, recKeys, rec_spz);
        if (quadKeys.size() > 0) {
            // write quad info
            std::ofstream fs;
            fs.open(io::gOutputDirectory + "/WJ_BOX/SOLID_QUAD_KEY_NR.rank"
                    + mpi::strRank() + ".txt" );
            for (int iq = 0; iq < quadKeys.size(); iq++) {
                fs << quadKeys[iq] << " " << quadNrs[iq] << "\n";
            }
            fs.close();
            
            // write receivers
            fs.open(io::gOutputDirectory + "/WJ_BOX/SOLID_REC_KEY_SPZ.rank"
                    + mpi::strRank() + ".txt" );
            fs.precision(18);
            for (int iq = 0; iq < quadKeys.size(); iq++) {
                for (int irec = 0; irec < recKeys[iq].size(); irec++) {
                    fs << recKeys[iq][irec] << " "
                    << rec_spz[iq].row(irec) << "\n";
                }
            }
            fs.close();
        }
        
        // fluid
        wj->infoFluid(quadKeys, quadNrs, quadBoundaryPnts, recKeys, rec_spz);
        if (quadKeys.size() > 0) {
            // write quad info
            std::ofstream fs;
            fs.open(io::gOutputDirectory + "/WJ_BOX/FLUID_QUAD_KEY_NR.rank"
                    + mpi::strRank() + ".txt" );
            for (int iq = 0; iq < quadKeys.size(); iq++) {
                fs << quadKeys[iq] << " " << quadNrs[iq] << "\n";
            }
            fs.close();
            
            // write receivers
            fs.open(io::gOutputDirectory + "/WJ_BOX/FLUID_REC_KEY_SPZ.rank"
                    + mpi::strRank() + ".txt" );
            fs.precision(18);
            for (int iq = 0; iq < quadKeys.size(); iq++) {
                for (int irec = 0; irec < recKeys[iq].size(); irec++) {
                    fs << recKeys[iq][irec] << " "
                    << rec_spz[iq].row(irec) << "\n";
                }
            }
            fs.close();
        }
        
        // quit if info only
        mpi::barrier();
        bool infoOnly;
        wjfs >> prompt >> infoOnly;
        if (infoOnly) {
            exit(0);
        }
        
        // stf
        std::string fnameSTF_S, fnameSTF_F;
        wjfs >> prompt >> fnameSTF_S >> fnameSTF_F    ;
        bool aligned;
        wjfs >> prompt >> aligned;
        int bufferSize;
        wjfs >> prompt >> bufferSize;
        wj->initializeSTFs(io::gInputDirectory + "/" + fnameSTF_S,
			   io::gInputDirectory+"/"+fnameSTF_F,
                           aligned, bufferSize);
	if (mpi::root()){
	io::cout<<io::gInputDirectory+"/"+fnameSTF_S<<io::endl;
	io::cout<<io::gInputDirectory+"/"+fnameSTF_F<<io::endl;
	}
        domain->setWavefieldInjection(wj);
    }
    
    
    //////////////// source ////////////////
    bool hasSource = ifs.get<bool>();
    if (hasSource) {
        // element
        int etag = ifs.get<int>();
        auto srcEle = domain->getSolidElement(etag);
        // stf
        const std::vector<double> &times = ifs.get<std::vector<double>>();
        const eigen::RColX &singleSTF = ifs.get<eigen::RColX>();
        const eigen::CMatXN3 &singlePattern = ifs.get<eigen::CMatXN3>();
        auto stf = std::make_unique<const SolidForce::STF>
        (true, times, singleSTF, singlePattern);
        // element source
        std::unique_ptr<const ElementSource> force =
        std::make_unique<const SolidForce>(srcEle, stf);
        domain->addElementSource(force);
    }
    
    //////////////// time scheme ////////////////
    double t0 = ifs.get<double>();
    double dt = ifs.get<double>();
    int numTimeSteps = ifs.get<int>();
    double user_starting_time;
    pfs >> prompt >> user_starting_time;
    if (user_starting_time > t0) {
        // Compute n_start
        int n_start = std::ceil((user_starting_time - t0) / dt);

        // update
        t0 = t0 + n_start * dt;
        numTimeSteps = numTimeSteps - n_start;
    }
    timeScheme->setTime(t0, dt, numTimeSteps);
    
    //////////////// stations ////////////////
    // default values
    int numRecordSteps = ifs.get<int>();
    int sampleIntv = ifs.get<int>();
    int dumpIntv = ifs.get<int>();
    std::string comp = ifs.get<std::string>();
    channel::WavefieldCS wcs;
    if (comp == "SPZ") {
        wcs = channel::WavefieldCS::spz;
    } else if (comp == "RTZ") {
        wcs = channel::WavefieldCS::RTZ;
    } else {
        wcs = channel::WavefieldCS::ENZ;
    }
    
    // group
    int excludeUngrouped, ngroup;
    pfs >> prompt >> excludeUngrouped;
    pfs >> prompt >> ngroup;
    
    // group name
    std::vector<std::string> grpNames(ngroup + 2, "");
    pfs >> prompt;
    for (int ig = 0; ig < ngroup; ig++) {
        pfs >> grpNames[ig];
    }
    grpNames[ngroup] = "default_solid";
    grpNames[ngroup + 1] = "default_fluid";
    
    // group fluid
    std::vector<int> grpFluid(ngroup + 2, 0);
    pfs >> prompt;
    for (int ig = 0; ig < ngroup; ig++) {
        pfs >> grpFluid[ig];
    }
    grpFluid[ngroup] = 0;
    grpFluid[ngroup + 1] = 1;
    
    // group channel
    std::vector<std::vector<std::string>> grpChs(ngroup + 2,
                                                 std::vector<std::string>());
    pfs >> prompt;
    for (int ig = 0; ig < ngroup; ig++) {
        std::string chplus;
        pfs >> chplus;
        grpChs[ig] = bstring::split(chplus, "+");
    }
    grpChs[ngroup] = {"U"};
    grpChs[ngroup + 1] = {"U"};
    
    // group sample
    std::vector<int> grpSampleIntv(ngroup + 2, 0);
    std::vector<int> grpRecordSteps(ngroup + 2, 0);
    pfs >> prompt;
    for (int ig = 0; ig < ngroup; ig++) {
        pfs >> grpSampleIntv[ig];
        if (grpSampleIntv[ig] <= 0) {
            grpSampleIntv[ig] = sampleIntv;
        }
        grpRecordSteps[ig] = numTimeSteps / grpSampleIntv[ig];
        if (numTimeSteps % grpSampleIntv[ig] > 0) {
            grpRecordSteps[ig]++;
        }
    }
    grpSampleIntv[ngroup] = sampleIntv;
    grpSampleIntv[ngroup + 1] = sampleIntv;
    grpRecordSteps[ngroup] = numRecordSteps;
    grpRecordSteps[ngroup + 1] = numRecordSteps;
    
    // dump
    std::vector<int> grpDumpIntv(ngroup + 2, 0);
    pfs >> prompt;
    for (int ig = 0; ig < ngroup; ig++) {
        pfs >> grpDumpIntv[ig];
        if (grpDumpIntv[ig] <= 0) {
            grpDumpIntv[ig] = dumpIntv;
        }
    }
    grpDumpIntv[ngroup] = dumpIntv;
    grpDumpIntv[ngroup + 1] = dumpIntv;
    
    // wcs
    std::vector<channel::WavefieldCS> grpWCS(ngroup + 2,
                                             channel::WavefieldCS::ENZ);
    pfs >> prompt;
    for (int ig = 0; ig < ngroup; ig++) {
        std::string comp;
        pfs >> comp;
        if (comp == "SPZ") {
            grpWCS[ig] = channel::WavefieldCS::spz;
        } else if (comp == "RTZ") {
            grpWCS[ig] = channel::WavefieldCS::RTZ;
        } else if (comp == "ENZ") {
            grpWCS[ig] = channel::WavefieldCS::ENZ;
        } else {
            grpWCS[ig] = wcs;
        }
    }
    grpWCS[ngroup] = wcs;
    grpWCS[ngroup + 1] = wcs;
    
    // groups
    std::vector<std::unique_ptr<StationGroup<StationSolid>>> solidGroups;
    std::vector<std::unique_ptr<StationGroup<StationFluid>>> fluidGroups;
    for (int ig = 0; ig < ngroup + 2; ig++) {
        // io
        std::unique_ptr<StationIO> stio = std::make_unique<StationIO_NetCDF>();
        // group
        std::vector<std::string> ch = grpChs[ig];
        if (grpFluid[ig]) {
            fluidGroups.
            push_back(std::make_unique<StationGroup<StationFluid>>
                      (grpNames[ig], grpRecordSteps[ig], grpSampleIntv[ig],
                       grpDumpIntv[ig], grpWCS[ig], ch, stio));
        } else {
            solidGroups.
            push_back(std::make_unique<StationGroup<StationSolid>>
                      (grpNames[ig], grpRecordSteps[ig], grpSampleIntv[ig],
                       grpDumpIntv[ig], grpWCS[ig], ch, stio));
        }
    }
    
    int nst = ifs.get<int>();
    for (int ist = 0; ist < nst; ist++) {
        // station
        std::string nwk = ifs.get<std::string>();
        std::string name = ifs.get<std::string>();
        std::string key = nwk + "." + name;
        double phi = ifs.get<double>();
        double theta = ifs.get<double>();
        double baz = ifs.get<double>();
        bool inFluid = ifs.get<bool>();
        if (inFluid) {
            std::unique_ptr<StationFluid> station =
            std::make_unique<StationFluid>(key, phi, theta, baz);
            // element + weights
            int etag = ifs.get<int>();
            auto stEle = domain->getFluidElement(etag);
            eigen::DRowN weights = ifs.get<eigen::DRowN>();
            station->setElement(stEle, weights);
            // find group by network
            for (int ig = 0; ig < fluidGroups.size(); ig++) {
                if (nwk == fluidGroups[ig]->getGroupName() ||
                    (fluidGroups[ig]->getGroupName() == "default_fluid"
                     && !excludeUngrouped)) {
                    fluidGroups[ig]->addStation(station);
                    break;
                }
            }
        } else {
            std::unique_ptr<StationSolid> station =
            std::make_unique<StationSolid>(key, phi, theta, baz);
            // element + weights
            int etag = ifs.get<int>();
            auto stEle = domain->getSolidElement(etag);
            eigen::DRowN weights = ifs.get<eigen::DRowN>();
            station->setElement(stEle, weights);
            // find group by network
            for (int ig = 0; ig < solidGroups.size(); ig++) {
                if (nwk == solidGroups[ig]->getGroupName() ||
                    (solidGroups[ig]->getGroupName() == "default_solid"
                     && !excludeUngrouped)) {
                    solidGroups[ig]->addStation(station);
                    break;
                }
            }
        }
    }
    
    // add to domain
    for (int ig = 0; ig < solidGroups.size(); ig++) {
        domain->addStationGroupInSolid(solidGroups[ig]);
    }
    for (int ig = 0; ig < fluidGroups.size(); ig++) {
        domain->addStationGroupInFluid(fluidGroups[ig]);
    }
    
    // close bin file
    ifs.close();
    
    //////////////// static ////////////////
    GradientQuadrature<numerical::Real>::
    setGMat(se_constants::gGMatrixGLL, se_constants::gGMatrixGLJ);
    
    // fftw
    fftw::createPlans(5.);
    if (io::gVerbose == io::VerboseLevel::detailed) {
        io::cout << fftw::verbose();
    }
    
    //////////////// solve ////////////////
    // solve
    timeScheme->solve();
    
    // fftw
    fftw::clearPlans();
}
