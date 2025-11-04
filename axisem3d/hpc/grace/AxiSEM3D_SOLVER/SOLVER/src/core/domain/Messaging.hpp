//
//  Messaging.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 4/29/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  mpi messaging

#ifndef Messaging_hpp
#define Messaging_hpp

#include "MessageRank.hpp"

class Messaging {
public:
    // add a rank to communicate with
    void addRankComm(std::unique_ptr<MessageRank> &msgRank) {
        mMessageRanks.push_back(std::move(msgRank));
        mRequestsSend.push_back(MPI_REQUEST_NULL);
        mRequestsRecv.push_back(MPI_REQUEST_NULL);
    }
    
    // phase 1: gather, send, recv
    void commGatherSendRecv() {
        for (int icomm = 0; icomm < mMessageRanks.size(); icomm++) {
            // gather
            mMessageRanks[icomm]->gatherFromPoints();
            // send
            mMessageRanks[icomm]->sendBuffer(mRequestsSend[icomm]);
            // recv
            mMessageRanks[icomm]->recvBuffer(mRequestsRecv[icomm]);
        }
    }
    
    // phase 2: wait, scatter
    void commWaitScatter() {
        // wait for recv
        mpi::wait_all((int)mRequestsRecv.size(), mRequestsRecv.data());
        // scatter
        for (int icomm = 0; icomm < mMessageRanks.size(); icomm++) {
            mMessageRanks[icomm]->scatterToPoints();
        }
        // wait for send
        mpi::wait_all((int)mRequestsSend.size(), mRequestsSend.data());
    }
    
    // check if a point exists in a domain of smaller rank
    bool pointInSmallerRank(const std::shared_ptr<const Point> &target) const {
        for (int icomm = 0; icomm < mMessageRanks.size(); icomm++) {
            if (mMessageRanks[icomm]->getRankOther() < mpi::rank()
                && mMessageRanks[icomm]->contains(target)) {
                return true;
            }
        }
        return false;
    }
    
    // get size
    int getNumRankComm() const {
        return (int)mMessageRanks.size();
    }
    
private:
    // message on ranks
    std::vector<std::unique_ptr<MessageRank>> mMessageRanks;
    
    // requests
    std::vector<MPI_Request> mRequestsSend;
    std::vector<MPI_Request> mRequestsRecv;
    
    
    ////////////////////////////////////////
    //////////////// static ////////////////
    ////////////////////////////////////////
    
public:
    // stream-based creator
    static std::unique_ptr<Messaging>
    createFromStream(sbs::ifstream &ifs, const Domain &domain) {
        std::unique_ptr<Messaging> msg = std::make_unique<Messaging>();
        int ncom = ifs.get<int>();
        for (int icom = 0; icom < ncom; icom++) {
            // create and add message on a single rank
            auto msgRank = std::make_unique<MessageRank>(ifs, domain);
            msg->addRankComm(msgRank);
        }
        return msg;
    }
};

#endif /* Messaging_hpp */
