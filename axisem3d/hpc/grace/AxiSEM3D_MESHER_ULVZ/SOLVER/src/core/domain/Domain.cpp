// Domain.cpp
// created by Kuangdai on 8-Apr-2016 
// computational domain

#include "Domain.h"
#include "Point.h"
#include "SolidPoint.h"
#include "FluidPoint.h"
#include "Element.h"
#include "SolidElement.h"
#include "FluidElement.h"
#include "SolidFluidPoint.h"
#include "SourceTerm.h"
#include "SourceTimeFunction.h"
#include "PointwiseRecorder.h"
#include "SurfaceRecorder.h"
#include "XMPI.h"
#include "NuWisdom.h"
#include "MultilevelTimer.h"

Domain::Domain() {
    #ifdef _MEASURE_TIMELOOP
        mTimerElemts = new MyBoostTimer();
        mTimerPoints = new MyBoostTimer();
        mTimerAssemb = new MyBoostTimer();
        mTimerAsWait = new MyBoostTimer();
        mTimerOthers = new MyBoostTimer();
    #endif
}

Domain::~Domain() {
    for (const auto &e: mPoints) {delete e;}
    for (const auto &e: mElements) {delete e;}
    for (const auto &e: mSourceTerms) {delete e;}
    if (mPointwiseRecorder) {delete mPointwiseRecorder;};
    if (mSurfaceRecorder) {delete mSurfaceRecorder;};
    if (mSTF) {delete mSTF;}
    if (mMsgInfo) {delete mMsgInfo;}
    if (mMsgBuffer) {delete mMsgBuffer;}
    if (mLearnPar) {delete mLearnPar;}
    #ifdef _MEASURE_TIMELOOP
        delete mTimerElemts;
        delete mTimerPoints;
        delete mTimerAssemb;
        delete mTimerAsWait;
        delete mTimerOthers;
    #endif
}

int Domain::addPoint(Point *point) {
    point->setDomainTag(mPoints.size());
    mPoints.push_back(point);
    return point->getDomainTag();
}


int Domain::addSolidPoint(SolidPoint *point) {
    point->mSolidTag=mSolidPoints.size();
    mSolidPoints.push_back(point);
    return point->mSolidTag;
}

int Domain::addFluidPoint(FluidPoint *point) {
    point->mFluidTag=mFluidPoints.size();
    mFluidPoints.push_back(point);
    return point->mFluidTag;
}

int Domain::addElement(Element *elem) {
    elem->setDomainTag(mElements.size());
    mElements.push_back(elem);
    return elem->getDomainTag();
}

int Domain::addSolidElement(SolidElement *elem) {
    elem->mSolidTag=(mSolidElements.size());
    mSolidElements.push_back(elem);
    return elem->mSolidTag;
}
int Domain::addFluidElement(FluidElement *elem) {
    elem->mFluidTag=(mFluidElements.size());
    mFluidElements.push_back(elem);
    return elem->mFluidTag;
}

void Domain::test() const {
    for (const auto &point: mPoints) {
        point->test();
    }
    for (const auto &elem: mElements) {
        elem->test();
    }
}

void Domain::resetZero() const {
    for (const auto &point: mPoints) {
        point->resetZero();
    }
    for (const auto &elem: mElements) {
        elem->resetZero();
    }
}

void Domain::initDisplTinyRandom() const {
    for (const auto &point: mPoints) {
        point->randomDispl((Real)1e-30, point->getDomainTag(), 2);
    }
}

void Domain::computeStiff() const {
    #ifdef _MEASURE_TIMELOOP
        mTimerElemts->resume();
    #endif
    
    for (const auto &elem: mElements) {
        elem->computeStiff();
    }
    
    #ifdef _MEASURE_TIMELOOP
        mTimerElemts->stop();
    #endif
}

void Domain::applySource(int tstep) const {
    #ifdef _MEASURE_TIMELOOP
        mTimerElemts->resume();
    #endif
    
    Real stf = mSTF->getFactor(tstep);
    for (const auto &source: mSourceTerms) {
        source->apply(stf);
    }
    
    #ifdef _MEASURE_TIMELOOP
        mTimerElemts->stop();
    #endif
}

void Domain::assembleStiff(int phase) const {
    #ifdef _MEASURE_TIMELOOP
        mTimerAssemb->resume();
    #endif
    
    if (phase <= 0) {
        // feed buffer
        for (int i = 0; i < mMsgInfo->mNProcComm; i++) {
            int row = 0;
            for (int j = 0; j < mMsgInfo->mNLocalPoints[i]; j++) {
                int pTag = mMsgInfo->mILocalPoints[i][j];
                mPoints[pTag]->feedBuffer(mMsgBuffer->mBufferSend[i], row);
            }
        }
        
        // send and recv
        for (int i = 0; i < mMsgInfo->mNProcComm; i++) {
            XMPI::isendComplex(mMsgInfo->mIProcComm[i], mMsgBuffer->mBufferSend[i], mMsgInfo->mReqSend[i]);
            XMPI::irecvComplex(mMsgInfo->mIProcComm[i], mMsgBuffer->mBufferRecv[i], mMsgInfo->mReqRecv[i]);
        }
    }
    
    if (phase >= 0) {
        // extract buffer 
        #ifdef _MEASURE_TIMELOOP
            mTimerAsWait->resume();
        #endif
        XMPI::wait_all(mMsgInfo->mReqRecv.size(), mMsgInfo->mReqRecv.data());
        #ifdef _MEASURE_TIMELOOP
            mTimerAsWait->stop();
        #endif
        
        for (int i = 0; i < mMsgInfo->mNProcComm; i++) {
            int row = 0;
            for (int j = 0; j < mMsgInfo->mNLocalPoints[i]; j++) {
                int pTag = mMsgInfo->mILocalPoints[i][j];
                mPoints[pTag]->extractBuffer(mMsgBuffer->mBufferRecv[i], row);
            }
        }
        
        #ifdef _MEASURE_TIMELOOP
            mTimerAsWait->resume();
        #endif
        XMPI::wait_all(mMsgInfo->mReqSend.size(), mMsgInfo->mReqSend.data());
        #ifdef _MEASURE_TIMELOOP
            mTimerAsWait->stop();
        #endif
    }
    
    #ifdef _MEASURE_TIMELOOP
        mTimerAssemb->stop();
    #endif
}

void Domain::updateNewmark(double dt) const {
    #ifdef _MEASURE_TIMELOOP
        mTimerPoints->resume();
    #endif
    
    for (const auto &point: mPoints) {
        point->updateNewmark(dt);
    }
    
    #ifdef _MEASURE_TIMELOOP
        mTimerPoints->stop();
    #endif
}

void Domain::coupleSolidFluid() const {
    #ifdef _MEASURE_TIMELOOP
        mTimerPoints->resume();
    #endif
    
    for (const auto &point: mSFPoints) {
        point->coupleSolidFluid();
    }
    
    #ifdef _MEASURE_TIMELOOP
        mTimerPoints->stop();
    #endif
}

void Domain::initializeRecorders() const {
    mPointwiseRecorder->initialize();
    if (mSurfaceRecorder) {
        mSurfaceRecorder->initialize();
    }
}

void Domain::finalizeRecorders() const {
    mPointwiseRecorder->finalize();
    if (mSurfaceRecorder) {
        mSurfaceRecorder->finalize();
    }
}

void Domain::record(int tstep, double t) const {
    #ifdef _MEASURE_TIMELOOP
        mTimerOthers->resume();
    #endif
    
    mPointwiseRecorder->record(tstep, t);
    if (mSurfaceRecorder) {
        mSurfaceRecorder->record(tstep, t);
    }
    
    #ifdef _MEASURE_TIMELOOP
        mTimerOthers->stop();
    #endif
}

void Domain::dumpLeft() const {
    #ifdef _MEASURE_TIMELOOP
        mTimerOthers->resume();
    #endif
    
    mPointwiseRecorder->dumpToFile();
    if (mSurfaceRecorder) {
        mSurfaceRecorder->dumpToFile();
    }
    
    #ifdef _MEASURE_TIMELOOP
        mTimerOthers->stop();
    #endif
}

void Domain::checkStability(double dt, int tstep, double t) const {
    #ifdef _MEASURE_TIMELOOP
        mTimerOthers->resume();
    #endif
    
    bool unstable = false;
    Point *unstable_point = 0;
    
    for (const auto &point: mPoints) {
        if (!point->stable()) {
            unstable = true;
            unstable_point = point;
            break;
        }
    }
    
    #ifdef _MEASURE_TIMELOOP
        mTimerOthers->stop();
    #endif
    
    if (unstable) {
        const RDCol2 &sz = unstable_point->getCoords() / 1e3;
        double r = sz.norm();
        double theta = (r < tinyDouble) ? 0. : acos(sz(1) / r);
        XMPI::cout.setp(XMPI::rank());
        XMPI::cout << "\n*****************************************" << XMPI::endl;
        XMPI::cout << "  SIMULATION BLEW UP! AXISEM3D ABORTED!" << XMPI::endl << XMPI::endl;
        XMPI::cout << "  Where instability occured " << XMPI::endl;
        XMPI::cout << "    (s, z) / km    =   (" << sz(0) << ", " << sz(1) << ")" << XMPI::endl;
        XMPI::cout << "    Radius / km    =   " << r << XMPI::endl;
        XMPI::cout << "    Dist / degr    =   " << theta / degree << XMPI::endl;
        XMPI::cout << "  When instability occured " << XMPI::endl;
        XMPI::cout << "    Current time   =   " << t << XMPI::endl;
        XMPI::cout << "    Current step   =   " << tstep << XMPI::endl;
        XMPI::cout << "  DT = " << dt << " s. Try smaller ones." << XMPI::endl;
        XMPI::cout << "*****************************************\n" << XMPI::endl;
        throw std::runtime_error("Domain::checkStability || Simulation blew up.");
    }
}

#include <sstream>
#include <map>
#include <iomanip>
std::string Domain::verbose() const {
    // elements
    int nele = XMPI::sum((int)mElements.size());
    std::map<std::string, int> eles;
    for (const auto &elem: mElements) {
        std::string estr = elem->verbose();
        eles.insert(std::pair<std::string, int>(estr, 0));
        eles.at(estr) += 1;
    }
    
    std::vector<std::map<std::string, int>> all_eles;
    XMPI::gather(eles, all_eles, MPI_INT, false);
    if (XMPI::root()) {
        for (int iproc = 1; iproc < XMPI::nproc(); iproc++) {
            for (auto it = all_eles[iproc].begin(); it != all_eles[iproc].end(); it++) {
                std::string estr = it->first;
                eles.insert(std::pair<std::string, int>(estr, 0));
                eles.at(estr) += it->second;
            }
        }
    }
    
    // points 
    int npoint = XMPI::sum((int)mPoints.size());
    std::map<std::string, int> points;
    for (const auto &point: mPoints) {
        std::string estr = point->verbose();
        points.insert(std::pair<std::string, int>(estr, 0));
        points.at(estr) += 1;
    }
    
    std::vector<std::map<std::string, int>> all_points;
    XMPI::gather(points, all_points, MPI_INT, false);
    if (XMPI::root()) {
        for (int iproc = 1; iproc < XMPI::nproc(); iproc++) {
            for (auto it = all_points[iproc].begin(); it != all_points[iproc].end(); it++) {
                std::string estr = it->first;
                points.insert(std::pair<std::string, int>(estr, 0));
                points.at(estr) += it->second;
            }
        }
    }
    
    std::stringstream ss;
    ss << "\n=================== Computational Domain ===================" << std::endl;
    ss << "  Elements__________________________________________________" << std::endl;
    int width = 12;
    for (auto it = eles.begin(); it != eles.end(); it++) {
        width = std::max(width, (int)(it->first.length()));
    }
    ss << "    " << std::setw(width) << std::left << "TOTAL NUMBER" << "   =   " << nele << std::endl;
    for (auto it = eles.begin(); it != eles.end(); it++) {
        ss << "    " << std::setw(width) << std::left << it->first << "   =   " << it->second << std::endl;
    }

    ss << "  GLL Points________________________________________________" << std::endl;
    width = 12;
    for (auto it = points.begin(); it != points.end(); it++) {
        width = std::max(width, (int)(it->first.length()));
    }
    ss << "    " << std::setw(width) << std::left << "TOTAL NUMBER" << "   =   " << npoint << std::endl;
    for (auto it = points.begin(); it != points.end(); it++) {
        ss << "    " << std::setw(width) << std::left << it->first << "   =   " << it->second << std::endl;
    }
    
    ss << "=================== Computational Domain ===================\n" << std::endl;
    return ss.str();
}

std::string Domain::reportCost() const {
    #ifdef _MEASURE_TIMELOOP
        const double sec2h = 1. / 3600.;
        std::stringstream ss;
        double costE = mTimerElemts->elapsed() * sec2h;
        double costP = mTimerPoints->elapsed() * sec2h;
        double costA = mTimerAssemb->elapsed() * sec2h;
        double costW = mTimerAsWait->elapsed() * sec2h;
        double costO = mTimerOthers->elapsed() * sec2h;
        double costT = costE + costP + costA + costO;
        ss.precision(4);
        ss << std::setw(10) << std::left << XMPI::rank() << "   ";
        ss << std::setw(10) << std::left << costT << "  ";
        ss << std::setw(10) << std::left << costE << "      ";
        ss << std::setw(10) << std::left << costP << "    ";
        ss << std::setw(10) << std::left << costA - costW << "      ";
        ss << std::setw(10) << std::left << costW << "  ";
        ss << std::setw(10) << std::left << costO << std::endl;
        std::vector<std::string> all_s;
        XMPI::gather(ss.str(), all_s, false);
        
        std::string s = "";
        if (XMPI::root()) {
            s += "\n-------------------------------------- MPI COST MEASUREMENTS --------------------------------------\n";
            s += "PROCESSOR    WALLTIME    ELEMENT-WISE    POINT_WISE    MPI_ASSEMBLE    MPI_WAIT    MISCELLANEOUS\n";
            for (int iproc = 0; iproc < XMPI::nproc(); iproc++) {
                s += all_s[iproc]; 
            }
            s += "---------------------------------------------------------------------------------------------------\n\n\n";
        }
        return s;
    #else 
        return "";
    #endif
}

void Domain::learnWisdom(int tstep) const {
    if (!mLearnPar->mInvoked) {
        return;
    }
    
    #ifdef _MEASURE_TIMELOOP
        mTimerOthers->resume();
    #endif
    
    if (tstep % mLearnPar->mInterval == 0) {
        for (const auto &point: mPoints) {
            point->learnWisdom(mLearnPar->mCutoff);
        }
    }
    
    #ifdef _MEASURE_TIMELOOP
        mTimerOthers->stop();
    #endif
}

void Domain::dumpWisdom() const {
    if (!mLearnPar->mInvoked) {
        return;
    }
    
    #ifdef _MEASURE_TIMELOOP
        mTimerOthers->resume();
    #endif
    
    std::vector<double> buffer;
    for (const auto &point: mPoints) {
        if (!pointInPreviousRank(point->getDomainTag())) {
            buffer.push_back(point->getCoords()(0));
            buffer.push_back(point->getCoords()(1));
            buffer.push_back(point->getNuWisdom());
            buffer.push_back(point->getNu());
        }
    }
    
    std::vector<std::vector<double>> all_buffer;
    XMPI::gather(buffer, all_buffer, MPI_DOUBLE, false);
    if (XMPI::root()) {
        for (int iproc = 1; iproc < XMPI::nproc(); iproc++) {
            buffer.insert(buffer.end(), all_buffer[iproc].begin(), all_buffer[iproc].end());
        }
        NuWisdom wis;
        for (int i = 0; i < buffer.size() / 4; i++) {
            wis.insert(buffer[i * 4], buffer[i * 4 + 1], 
                round(buffer[i * 4 + 2]), round(buffer[i * 4 + 3]));
        }
        wis.writeToFile(mLearnPar->mFileName);
    }
    
    #ifdef _MEASURE_TIMELOOP
        mTimerOthers->stop();
    #endif
}

bool Domain::pointInPreviousRank(int myPointTag) const {
    for (int i = 0; i < mMsgInfo->mNProcComm; i++) {
        for (int j = 0; j < mMsgInfo->mNLocalPoints[i]; j++) {
            int pTag = mMsgInfo->mILocalPoints[i][j];
            if (pTag == myPointTag) {
                int otherRank = mMsgInfo->mIProcComm[i];
                if (otherRank < XMPI::rank()) {
                    return true;
                } else {
                    break;
                }
            }
        }
    }
    return false;
}

#include <iostream>

void Domain::dumpToStream(sbs::ofstream &ofs) const {
    // points
    int npntsolid = mSolidPoints.size();
    ofs << npntsolid;
    for (int ipnt = 0; ipnt < npntsolid; ipnt++) {
        mSolidPoints[ipnt]->dumpToStream(ofs);
    }
    
    int npntfluid = mFluidPoints.size();
    ofs << npntfluid;
    for (int ipnt = 0; ipnt < npntfluid; ipnt++) {
        mFluidPoints[ipnt]->dumpToStream(ofs);
    }
        
    int nelesolid = mSolidElements.size();
    ofs << nelesolid;
    for (int ie = 0; ie < nelesolid; ie++) {
        mSolidElements[ie]->dumpToStream(ofs);
    }
    
    int nelefluid = mFluidElements.size();
    ofs << nelefluid;
    for (int ie = 0; ie < nelefluid; ie++) {
        mFluidElements[ie]->dumpToStream(ofs);
    }

    // solid-fluid
    ofs << (int)(mSFPoints.size());
    for (int isf = 0; isf < mSFPoints.size(); isf++) {
        mSFPoints[isf]->dumpToStream(ofs);
    }
    
    // abc
    ofs << 0 << 0;
    
    // axial
    std::vector<int> atags;
    for (int ipnt = 0; ipnt < npntsolid; ipnt++) {
        if(mSolidPoints[ipnt]->axial()) {
            atags.push_back(ipnt);
        }
    }
    ofs << atags;
    
    atags.clear();
    for (int ipnt = 0; ipnt < npntfluid; ipnt++) {
        if(mFluidPoints[ipnt]->axial()) {
            atags.push_back(ipnt);
        }
    }
    ofs << atags;
    
    
    //  fuild surf
    ofs << 0;
    
    // mpi
    ofs<<mMsgInfo->mNProcComm;
    
    
    
    for (int i = 0; i < mMsgInfo->mNProcComm; i++) {
        int otherRank = mMsgInfo->mIProcComm[i];
        ofs<<otherRank;
        IMatXX table(mMsgInfo->mNLocalPoints[i], 3);
        for (int j = 0; j < mMsgInfo->mNLocalPoints[i]; j++) {
            int meshTag = mMsgInfo->mILocalPoints[i][j];
            table(j,0) = meshTag;
            table(j,1) = mMapMeshSolidFluid(meshTag, 1);
            table(j,2) = mMapMeshSolidFluid(meshTag, 2);
        }
        ofs<<table;
    }
    
    for (int irank = 0 ; irank<XMPI::nproc(); irank++) {
        if (XMPI::rank()==irank) {
            for (int i = 0; i < mMsgInfo->mNProcComm; i++) {
                int otherRank = mMsgInfo->mIProcComm[i];
                std::cout<<XMPI::rank()<<" => "<<otherRank <<" np:"
                << mMsgInfo->mNLocalPoints[i]<<std::endl;
            }
        }
        XMPI::barrier();
    }
    XMPI::barrier();
    
    
    ///////////////////////////
    // source
    bool hasSource = mSourceTerms.size()>0;
    ofs <<hasSource;
    if (hasSource) {
        
        ofs << mSourceTerms[0]->mElement->mSolidTag;
        auto time =std::vector<double>(mSTF->getSize(), 0);
        for (int i=0;i<time.size();i++) {
            time[i] = mSTF->mDeltaT * i - mSTF->mShift;
        }
        ofs<< time;
        
        RColX b = Eigen::Map<RColX>(mSTF->mSTF.data(), mSTF->getSize());
        ofs<< b;
        
        CMatXN3 force=CMatXN3::Zero(3, 75);
        for (int alpha = 0; alpha<=2;alpha++) {
            for (int idim=0;idim<3;idim++) {
                for (int ipol = 0;ipol<nPntEdge;ipol++) {
                    for (int jpol = 0;jpol<nPntEdge;jpol++) {
                        force(alpha, idim*25+ipol*5+jpol) =
                        mSourceTerms[0]->mForce[ipol*5+jpol](alpha, idim);
                    }
                }
            }
        }
        ofs << force;
    }
    XMPI::barrier();
    
    ofs<<-mSTF->mShift<<mSTF->mDeltaT<<mSTF->getSize();
    
    // stations
    int nstep = mPointwiseRecorder->mTotalRecordSteps;
    int sampleIntv= mPointwiseRecorder->mRecordInterval;
    int dumpIntv = mPointwiseRecorder->mBufferSize;
    ofs<<nstep<<sampleIntv<<dumpIntv;
    ofs<<mPointwiseRecorder->mComponents;
    int nst = mPointwiseRecorder->mPointwiseInfo.size();
    ofs<<nst;
    for(int ist=0;ist<nst;ist++) {
        ofs<<mPointwiseRecorder->mPointwiseInfo[ist].mNetwork;
        ofs<<mPointwiseRecorder->mPointwiseInfo[ist].mName;
        ofs << mPointwiseRecorder->mPointwiseInfo[ist].mPhi;
        ofs << mPointwiseRecorder->mPointwiseInfo[ist].mTheta;
        ofs << mPointwiseRecorder->mPointwiseInfo[ist].mBAz;
        
        if (mPointwiseRecorder->mPointwiseInfo[ist].mElement->verbose()
        .find("FluidElement") != std::string::npos) {
            ofs << true;
            ofs << mPointwiseRecorder->mPointwiseInfo[ist].mElement->mFluidTag;
        } else {
            ofs << false;
            ofs << mPointwiseRecorder->mPointwiseInfo[ist].mElement->mSolidTag;
        }
        RDRowN w;
        for (int ipol = 0;ipol<nPntEdge;ipol++) {
            for (int jpol = 0;jpol<nPntEdge;jpol++) {
                w(ipol*5+jpol)=mPointwiseRecorder->mPointwiseInfo[ist].mWeights(ipol,jpol);
            }
        }
        ofs << w;
    }
    
    
    ofs.close();
    XMPI::barrier();
    std::cout << "RANK " << XMPI::rank() << " FINISHES\n";
    exit(0);
}

