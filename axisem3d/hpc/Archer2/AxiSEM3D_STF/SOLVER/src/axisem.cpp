// axisem.cpp
// created by Kuangdai on 26-Mar-2016 
// main of AxiSEM3D

#include "axisem.h"
#include "MultilevelTimer.h"

#include "XMPI.h"
#include "eigenc.h"
#include "eigenp.h"
#include "sbs.hpp"

#include <set>
#include <iostream>
#include "Parameters.h"
#include "Geodesy.h"
#include <map>
#include "NetCDF_Reader.h"
#include "NetCDF_Writer.h"

#include <unsupported/Eigen/CXX11/Tensor>
#include "SolverFFTW_1_nonstatic.h"

#ifdef _NO_OMP
int omp_get_max_threads() {
    return 1;
}
int omp_get_num_threads() {
    return 1;
}
int omp_get_thread_num() {
    return 0;
}
#else
#include <omp.h>
#endif

RDCol3 getPoint(const RDCol3 &spz, double lat_src, double lon_src) {
    RDCol3 rtpS;
    rtpS(0) = sqrt(spz(0) * spz(0) + spz(2) * spz(2));
    rtpS(1) = acos(spz(2) / rtpS(0));
    rtpS(2) = spz(1);
    return Geodesy::rotateSrc2Glob(rtpS, lat_src, lon_src, 0.);
}

RDMat33 get_xspz_from_Global(const RDCol3 &rtpG, double slat, double slon) {
    RDCol3 rtp = Geodesy::rotateGlob2Src(rtpG, slat, slon, 0.);
    RDCol3 spz_me, spz_ds, spz_dz;
    spz_me(0) = rtp(0) * sin(rtp(1));
    spz_me(1) = rtp(2);
    spz_me(2) = rtp(0) * cos(rtp(1));
    spz_ds = spz_me;
    spz_dz = spz_me;
    spz_ds(0) += 100;
    spz_dz(2) += 100;
    const RDCol3 &rtp_ds = getPoint(spz_ds, slat, slon);
    const RDCol3 &rtp_dz = getPoint(spz_dz, slat, slon);
    const RDCol3 &xyz_ds = Geodesy::toCartesian(rtp_ds);
    const RDCol3 &xyz_dz = Geodesy::toCartesian(rtp_dz);
    const RDCol3 &xyz_me = Geodesy::toCartesian(rtpG);
    RDCol3 vs = xyz_ds - xyz_me;
    RDCol3 vz = xyz_dz - xyz_me;
    vs /= vs.norm();
    vz /= vz.norm();
    RDCol3 vp = vz.cross(vs);
    RDMat33 frame;
    frame.row(0) = vs.transpose();
    frame.row(1) = vp.transpose();
    frame.row(2) = vz.transpose();
    return frame;
}

RDMat33 get_xspz_from_Src(const RDCol3 &rtp, const RDCol3 &rtpG, double slat, double slon) {
    RDCol3 spz_me, spz_ds, spz_dz;
    spz_me(0) = rtp(0) * sin(rtp(1));
    spz_me(1) = rtp(2);
    spz_me(2) = rtp(0) * cos(rtp(1));
    spz_ds = spz_me;
    spz_dz = spz_me;
    spz_ds(0) += 100;
    spz_dz(2) += 100;
    const RDCol3 &rtp_ds = getPoint(spz_ds, slat, slon);
    const RDCol3 &rtp_dz = getPoint(spz_dz, slat, slon);
    const RDCol3 &xyz_ds = Geodesy::toCartesian(rtp_ds);
    const RDCol3 &xyz_dz = Geodesy::toCartesian(rtp_dz);
    const RDCol3 &xyz_me = Geodesy::toCartesian(rtpG);
    RDCol3 vs = xyz_ds - xyz_me;
    RDCol3 vz = xyz_dz - xyz_me;
    vs /= vs.norm();
    vz /= vz.norm();
    RDCol3 vp = vz.cross(vs);
    RDMat33 frame;
    frame.row(0) = vs.transpose();
    frame.row(1) = vp.transpose();
    frame.row(2) = vz.transpose();
    return frame;
}
                    

int axisem_main(int argc, char *argv[]) {
    
    try {
        
        // variable sets
        PreloopVariables pl;
        SolverVariables sv;
        
        // initialize mpi
        XMPI::initialize(argc, argv);
        
        //////////////////////////////////////////
        //////////////////////////////////////////
        //////////////////////////////////////////
        //////////////////////////////////////////
        //////////////////////////////////////////
        //////////////////////////////////////////
        
        std::cout<<"num max threads: "<<omp_get_max_threads()<<std::endl;
        
        // GeoPoint
        double lat_src = boost::lexical_cast<double>(argv[1]);
        double lon_src = boost::lexical_cast<double>(argv[2]);

        double ulvz_lat = boost::lexical_cast<double>(argv[3]);
        double ulvz_lon = boost::lexical_cast<double>(argv[4]);
        
        // medium
        std::string medium = argv[5];
        int ndim = 1;
        if (medium == "SOLID")
            ndim = 3;
            
        // dir
        std::string dir_stations = argv[12];
        std::string dir_incident = argv[13];
        std::string dir_info = argv[14];
        std::string dir_stf = ".";

        // injection rec positions
        std::cout << "* injection rec positions" << std::endl;
        std::cout.flush();
        
        std::fstream fs(dir_info + "/output/WJ_BOX/"+medium+"_REC_KEY_SPZ.txt", std::fstream::in);
        std::vector<std::string> allLines;
        std::string line;
        while (getline(fs, line)) {
            allLines.push_back(line);
        }
        fs.close();
        
        int nrec_ulvz = allLines.size();
        std::cout << "number of recs in ulvz: " << nrec_ulvz << std::endl;
        std::cout.flush();
        
        std::map<std::string, RDCol3> rec_ulvz_loc_g;
        std::map<std::string, RDCol3> rec_ulvz_loc_s;
        for (int irec=0;irec<nrec_ulvz;irec++) {
            std::vector<std::string> strs = Parameters::splitString(allLines[irec], "\t ");
            double s = boost::lexical_cast<double>(strs[1]);
            double p = boost::lexical_cast<double>(strs[2]);
            double z = boost::lexical_cast<double>(strs[3]);
            RDCol3 rtpS;
            rtpS(0) = sqrt(s*s+z*z);
            rtpS(1) = acos(z/rtpS(0));
            rtpS(2) = p;
            rec_ulvz_loc_g.insert(std::pair<std::string,RDCol3>(strs[0],
                Geodesy::rotateSrc2Glob(rtpS, ulvz_lat, ulvz_lon, 0.)));
            rec_ulvz_loc_s.insert(std::pair<std::string,RDCol3>(strs[0],
                rtpS));
        }
        
        
        std::cout << "first rec in ulvz: " <<rec_ulvz_loc_g.begin()->first
        <<" " << rec_ulvz_loc_g.begin()->second.transpose() << std::endl;
        std::cout.flush();
        
        
        std::cout << "* reading station-rank file" << std::endl;
        std::cout.flush();
        fs.open(dir_incident+"/output/stations/INCIDENT_"+medium+"/rank_station.info");
        allLines.clear();
        while (getline(fs, line)) {
            allLines.push_back(line);
        }
        fs.close();

        std::set<int> all_ranks;
        std::map<std::string, int> st_ranks, st_orders;
        for (int ist=1;ist<allLines.size();ist++) {
            std::vector<std::string> strs = Parameters::splitString(allLines[ist], "\t ");
            st_ranks.insert(std::pair<std::string,int>(strs[1], boost::lexical_cast<int>(strs[0])));
            st_orders.insert(std::pair<std::string,int>(strs[1], boost::lexical_cast<int>(strs[2])));
            all_ranks.insert(boost::lexical_cast<int>(strs[0]));
        }
        std::cout << "first station rank: "<<st_ranks.begin()->first<<" "<<st_ranks.begin()->second<<std::endl;
        std::cout << "first station order: "<<st_orders.begin()->first<<" "<<st_orders.begin()->second<<std::endl;
        std::cout.flush();
        
        
        
        // open nc
        std::cout << "* opening nc" << std::endl;
        std::cout.flush();
        std::map<int, NetCDF_Reader *> nc_files;
        for (int irank : all_ranks) {
            std::stringstream fname;
            fname << dir_incident<<
            "/output/stations/INCIDENT_"<<medium<<"/axisem3d_synthetics.nc.rank" << irank;
            NetCDF_Reader *reader=new NetCDF_Reader();
            try {
                reader->open(fname.str());
                nc_files.insert(std::pair<int, NetCDF_Reader *>(irank, reader));
                std::cout << "opened "<< fname.str()<< std::endl;
                std::cout.flush();
            } catch (...){
                delete reader;
            }
        }
         
        
        // read time
        std::cout<<"* reading time"<<std::endl;
        std::cout.flush();
        RDColX times;
        nc_files.begin()->second->read1D("time_points", times);
        //*******************************************************//
        //*******************************************************//
        int torigin = boost::lexical_cast<int>(argv[6]);
        int npts = boost::lexical_cast<int>(argv[7]);
        npts = std::min(npts, (int)(times.size() - torigin));
        //*******************************************************//
        //*******************************************************//
        std::cout<<"time steps = "<< npts <<std::endl;
        std::cout.flush();
        
        
        
        
        
        

        
        // interpolation
        // read grid
        std::cout<<"* interpolation grid"<<std::endl;
        std::cout.flush();
        fs.open(dir_stations+"/grid_dist.txt");
        allLines.clear();
        while (getline(fs, line)) {
            allLines.push_back(line);
        }
        fs.close();
        RDColX grid_dists(allLines.size());
        for (int i=0;i<allLines.size();i++) {
            grid_dists(i) = boost::lexical_cast<double>(allLines[i]);
        }
        
        fs.open(dir_stations+"/grid_depth_"+medium+".txt");
        allLines.clear();
        while (getline(fs, line)) {
            allLines.push_back(line);
        }
        fs.close();
        RDColX grid_depths(allLines.size());
        for (int i=0;i<allLines.size();i++) {
            grid_depths(i) = boost::lexical_cast<double>(allLines[i]);
        }
        
        int nazim_cardinal = boost::lexical_cast<int>(argv[8]);
        int nazim_cardinal_ulvz = boost::lexical_cast<int>(argv[9]);
        int nazim_fourier = nazim_cardinal / 2 + 1;
        
        std::cout << "grid distance = " << grid_dists.size() << " "<< grid_dists(0)<<"~"<<grid_dists(grid_dists.size()-1)<<std::endl;
        std::cout << "grid depth = " << grid_depths.size() << " "<< grid_depths(0)<<"~"<<grid_depths(grid_depths.size()-1)<<std::endl;
        std::cout << "grid azimuth = " << nazim_cardinal << " "<< 0<<"~"<<360/nazim_cardinal*(nazim_cardinal-1)<<std::endl;
        
        // cut dist
        double distMin = boost::lexical_cast<double>(argv[10]);
        double distMax = boost::lexical_cast<double>(argv[11]);
        int dindex0 = std::lower_bound(grid_dists.begin(), grid_dists.end(), distMin) - grid_dists.begin() - 1;
        int dindex1 = std::upper_bound(grid_dists.begin(), grid_dists.end(), distMax) - grid_dists.begin();
        dindex0 = std::max(dindex0, 0);
        dindex1 = std::min(dindex1, (int)grid_dists.size() - 1);
        int dcount = dindex1 - dindex0 + 1;
        RDColX grid_dists_cut(dcount);
        for (int i =dindex0; i<=dindex1;i++) {
            grid_dists_cut[i-dindex0] = grid_dists[i];
        }
        grid_dists = grid_dists_cut;
        std::cout << "grid distance cut = " << grid_dists.size() << " "<< grid_dists(0)<<"~"<<grid_dists(grid_dists.size()-1)<<std::endl;

        
        
        
        
        
        
        
        
        
        // grid data
        std::cout<<"* interpolation data on grid"<<std::endl;
        std::cout.flush();
        Eigen::Tensor<Complex, 5, Eigen::RowMajor>
        u_cardinal_EVT(grid_dists.size(), grid_depths.size(),
        nazim_fourier, npts, ndim);
        u_cardinal_EVT.setConstant(std::numeric_limits<Complex>::quiet_NaN());
        
        std::vector<std::pair<int,int>> failed_idist_idepth;
        
        for (auto it_nc=nc_files.begin(); it_nc!=nc_files.end(); it_nc++) {
            // read data on this rank file
            int rank = it_nc->first;
            std::cout<<"** doing rank = "<< rank << std::endl;
            std::cout.flush();
            std::vector<Real> temp;
            std::vector<size_t> dims;
            it_nc->second->readMetaData("data", temp, dims);
            std::cout<<"dims="<<dims[0]<<"*"<<dims[1]<<"*"<<dims[2]<<"="<<temp.size()<<std::endl;
            std::cout.flush();
            Eigen::Tensor<Real, 3, Eigen::RowMajor> nc_data(npts, ndim, (int)dims[2]);
            for (int i = 0; i < npts; i++) {
                for (int j = 0; j < ndim; j++) {
                    for (int k = 0; k < dims[2]; k++) {
                        int pos = k + j * dims[2] + (i+torigin) * ndim * dims[2];
                        nc_data(i, j, k) = temp[pos];
                    }
                }
            }
            
            temp.clear();
            temp.shrink_to_fit();
            std::cout<<"done reading" << std::endl;
            std::cout.flush();
            
            #pragma omp parallel 
            {
                
                // fft
                SolverFFTW_1_nonstatic fft;
                for (int i = 0; i < omp_get_num_threads(); i++) {
                    if (i == omp_get_thread_num()) {
                        fft.initialize(nazim_cardinal);
                    }
                    #pragma omp barrier 
                }
                
                
                std::stringstream ss;
                Eigen::array<int, 3> start1, start2;
                Eigen::array<int, 3> count1, count2;
                Eigen::array<int, 1> iarray;
                Eigen::array<int, 5> start5, count5;
                Eigen::array<int, 2> iarray2 = {npts, ndim};
                
                #pragma omp for
                for (int idist=0; idist<grid_dists.size(); idist++) {
                    for (int idepth=0; idepth<grid_depths.size(); idepth++) {
                        ss.str("");
                        ss << "INCIDENT_"<<medium<<".ULVZ_"<<idist+dindex0 <<"_"<<idepth<<"_"<<0;
                        std::string st_key = ss.str();
                        int st_rank = st_ranks.at(st_key);
                        if (st_rank!=rank) {
                            continue;
                        }
                
                        bool success=true;
                        Eigen::Tensor<Real, 3, Eigen::RowMajor> u_card(
                        nazim_cardinal, npts, ndim);
                
                        for (int iazim =0;iazim<nazim_cardinal;iazim++) {
                            ss.str("");
                            ss << "INCIDENT_"<<medium<<".ULVZ_"<<idist+dindex0<<"_"<<idepth<<"_"<<iazim;
                            st_key = ss.str();
                            st_rank = st_ranks.at(st_key);
                            if (st_rank!=rank) {
                                success=false;
                                break;
                            }
                
                            int st_order = st_orders.at(st_key);
                            start1 = {iazim, 0, 0};
                            count1 = {1,npts,ndim};
                            start2 = {0,0,st_order};
                            count2 = {npts, ndim, 1};
                            u_card.slice(start1, count1).reshape(iarray2) = 
                            nc_data.slice(start2, count2).reshape(iarray2);
                        }
                
                        if (success) {
                            for (int itime=0; itime<npts; itime++) {
                                for (int idim=0; idim<ndim; idim++) {
                                    start1 = {0, itime, idim};
                                    count1 = {nazim_cardinal, 1, 1};
                                    iarray = {nazim_cardinal};
                                    Eigen::TensorMap<Eigen::Tensor<Real,1,Eigen::RowMajor>>
                                    (fft.getR2C_RMat().data(), nazim_cardinal)=
                                    u_card.slice(start1, count1).reshape(iarray);
                                    fft.computeR2C(nazim_cardinal);
                                    start5 = {idist, idepth, 0, itime, idim};
                                    count5 = {1, 1, nazim_fourier, 1, 1};
                                    iarray = {nazim_fourier};
                                    u_cardinal_EVT.slice(start5, count5).reshape(iarray)=
                                    Eigen::TensorMap<Eigen::Tensor<Complex,1,Eigen::RowMajor>>
                                    (fft.getR2C_CMat().data(), nazim_fourier);
                                }
                            }
                        } else {
                            failed_idist_idepth.push_back(std::pair<int, int>(idist, idepth));
                        }
                    }
                    if (omp_get_thread_num() == 0) {
                        std::cout << "Thread 0 done with dist " << idist << " out of " << grid_dists.size() << std::endl;
                        std::cout.flush();
                    }
                }
                fft.finalize();
            }
        }

        
        if (failed_idist_idepth.size() > 0) {
            std::cout << "Failed with " << failed_idist_idepth.size() 
            << " receivers on multiple ranks"<<std::endl;
            exit(0);
        } else {
            std::cout << "Done with interplation." <<std::endl;
        }
        
        
        
        
        
        // quad info
        std::cout<<"* quad info"<<std::endl;
        std::cout.flush();
        std::map<std::string,int> quad_nrs;
        fs.open(dir_info+"/output/WJ_BOX/"+medium+"_QUAD_KEY_NR.txt");
        while (getline(fs, line)) {
            std::vector<std::string> strs = Parameters::splitString(line, "\t ");
            quad_nrs.insert(std::pair<std::string,int>(strs[0],
                boost::lexical_cast<int>(strs[1])));
        }
        fs.close();
        std::cout<<"num quad = "<< quad_nrs.size() <<std::endl;
        
        
        
        // output data
        std::cout<<"* output nc"<<std::endl;
        std::cout.flush();
        NetCDF_Writer ncw;
        std::stringstream ssw;
        ssw<< dir_stf << "/output/injection_"<<medium<<".nc";
        ncw.open(ssw.str(),  true);
        std::vector<size_t> dims;
        dims.push_back(npts);
        ncw.defineVariable<double>("time_points", dims);
        RDColX times_use = times.block(torigin, 0, npts, 1);
        ncw.writeVariableWhole("time_points", times_use);

            
        // quad
        for (auto itq=quad_nrs.begin(); itq!=quad_nrs.end(); itq++) {
            std::string quad_key = itq->first;
            int quad_nr = itq->second;
            int ncol = 25 * ndim * (quad_nr / 2 + 1);
            std::vector<size_t> dims;
            dims.push_back(npts);
            dims.push_back(ncol);
            ncw.defineVariable<Real>(quad_key + "_RE", dims);
            ncw.defineVariable<Real>(quad_key + "_IM", dims);
        }
        
        
        std::cout<<"* quad loop"<<std::endl;
        std::cout.flush();
        
        
        
        for (auto itq=quad_nrs.begin(); itq!=quad_nrs.end(); itq++) {
            
            std::string quad_key = itq->first;
            int quad_nr = itq->second;
            
        
            // real data
            Eigen::Tensor<Real, 4, Eigen::RowMajor> u_spz_BOX_cardinal(
                npts, quad_nr, ndim, 25);
            
            u_spz_BOX_cardinal.setConstant(std::numeric_limits<Real>::quiet_NaN());
                
            //####### nr loop #######
            #pragma omp parallel 
            {
                Eigen::array<int, 3> iarray3 = {nazim_fourier, npts, ndim};
                Eigen::array<int, 4> start4;
                Eigen::array<int, 4> count4;
                Eigen::array<int, 3> start1, start2;
                Eigen::array<int, 3> count1, count2;
                Eigen::array<int, 1> iarray;
                Eigen::array<int, 5> start5, count5;
                Eigen::array<int, 2> iarray2 = {npts, ndim};
                
                #pragma omp for
                for (int ir =0;ir<quad_nr;ir++) {
                    //####### point loop #######
                    for (int ipnt=0;ipnt<25;ipnt++) {
                        // read data
                        std::stringstream ss;
                        ss<< quad_key << "__" << ipnt << '_' << ir;
                        std::string st_key = ss.str();
                        RDCol3 rtp = Geodesy::rotateGlob2Src(rec_ulvz_loc_g.at(st_key),
                            lat_src, lon_src, 0.);
                        
                        
                        // for rotation
                        const RDMat33 &x_EVT = get_xspz_from_Global(rec_ulvz_loc_g.at(st_key), lat_src, lon_src);
                        const RDMat33 &x_BOX = get_xspz_from_Src(rec_ulvz_loc_s.at(st_key),
                                                                 rec_ulvz_loc_g.at(st_key), ulvz_lat, ulvz_lon);
                        
                        
                        // dist
                        double dist = rtp[1];
                        int idist1 = std::upper_bound(grid_dists.begin(), grid_dists.end(), dist)-
                        grid_dists.begin();
                        if (idist1 == 0) {
                            idist1 = 1;
                        }
                        if (idist1 == grid_dists.size()) {
                            idist1 = grid_dists.size() - 1;
                        }
                        int idist0 = idist1 - 1;
                        double fdist1 = (dist - grid_dists[idist0]) / (grid_dists[idist1] - grid_dists[idist0]);
                        double fdist0 = 1. - fdist1;
                    
                        // depth
                        double depth = (6371e3-rec_ulvz_loc_g[st_key](0))/1e3;
                        int idep1 = std::upper_bound(grid_depths.begin(), grid_depths.end(), depth)-
                        grid_depths.begin();
                        if (idep1 == 0) {
                            idep1 = 1;
                        }
                        if (idep1 == grid_depths.size()) {
                            idep1 = grid_depths.size() - 1;
                        }
                        int idep0 = idep1 - 1;
                        double fdep1 = (depth - grid_depths[idep0]) / (grid_depths[idep1] - grid_depths[idep0]) ;
                        double fdep0 = 1. - fdep1;
                        
                        // inplane interp
                        Eigen::Tensor<Complex, 3, Eigen::RowMajor> u_card(
                            nazim_fourier, npts, ndim);
                        u_card.setZero();
                        count5={1,1,nazim_fourier,npts,ndim};
                        start5={idist0, idep0, 0,0,0};
                        u_card += fdist0 * fdep0 * u_cardinal_EVT.slice(start5,count5).reshape(iarray3);
                        start5={idist1, idep0, 0,0,0};
                        u_card += fdist1 * fdep0 * u_cardinal_EVT.slice(start5,count5).reshape(iarray3);
                        start5={idist0, idep1, 0,0,0};
                        u_card += fdist0 * fdep1 * u_cardinal_EVT.slice(start5,count5).reshape(iarray3);
                        start5={idist1, idep1, 0,0,0};
                        u_card += fdist1 * fdep1 * u_cardinal_EVT.slice(start5,count5).reshape(iarray3);
                        
                        // Fourier
                        double phi = rtp[2];
                        Eigen::Tensor<Real, 2, Eigen::RowMajor> u_spz_EVT(npts, ndim);
                        start1={0,0,0};
                        count1={1,npts,ndim};
                        u_spz_EVT.setZero();
                        u_spz_EVT += u_card.slice(start1, count1).reshape(iarray2).real();
                        for  (int alpha = 1; alpha < nazim_fourier; alpha++) {
                            start1={alpha,0,0};
                            u_spz_EVT += 2. *
                            ((Complex)exp(1.*alpha* phi * iid) *
                             u_card.slice(start1, count1).reshape(iarray2)).real();
                        }
                        
                        
                        if (ndim == 3) {
                            Eigen::Tensor<Real, 2, Eigen::RowMajor> u_spz_BOX(npts, ndim);
                            u_spz_BOX.setZero();
                            for (int j=0;j<3;j++) {
                                for (int i=0;i<3;i++) {
                                    double fact = x_EVT.row(i).dot(x_BOX.row(j));
                                    for (int itime=0;itime<npts;itime++) {
                                        u_spz_BOX(itime, j) += u_spz_EVT(itime, i) * fact;
                                    }
                                }
                            }                        
                            
                            start4={0,ir,0,ipnt};
                            count4={npts,1,ndim,1};
                            u_spz_BOX_cardinal.slice(start4,count4).reshape(iarray2) = u_spz_BOX;
                        } else {
                            start4={0,ir,0,ipnt};
                            count4={npts,1,ndim,1};
                            u_spz_BOX_cardinal.slice(start4,count4).reshape(iarray2) = u_spz_EVT;
                        }
                    }
                }
            }
            
            int quad_nc = quad_nr/2+1;
            Eigen::Tensor<Complex, 4, Eigen::RowMajor> u_spz_BOX_Fourier(
                npts, quad_nc, ndim, 25);
            u_spz_BOX_Fourier.setConstant(std::numeric_limits<Complex>::quiet_NaN());

            #pragma omp parallel 
            {
                Eigen::array<int, 3> iarray3 = {nazim_fourier, npts, ndim};
                Eigen::array<int, 4> start4;
                Eigen::array<int, 4> count4;
                Eigen::array<int, 3> start1, start2;
                Eigen::array<int, 3> count1, count2;
                Eigen::array<int, 1> iarray;
                Eigen::array<int, 5> start5, count5;
                Eigen::array<int, 2> iarray2 = {npts, ndim};
                
                // fft
                SolverFFTW_1_nonstatic fft;
                for (int i = 0; i < omp_get_num_threads(); i++) {
                    if (i == omp_get_thread_num()) {
                        fft.initialize(nazim_cardinal_ulvz);
                    }
                    #pragma omp barrier 
                }
                
                #pragma omp for
                for (int itime=0; itime<npts; itime++) {
                    for (int idim=0; idim<ndim; idim++) {
                        for (int ipnt=0; ipnt<25; ipnt++) {
                            start4 = {itime, 0, idim, ipnt};
                            count4 = {1, quad_nr, 1, 1};
                            iarray = {quad_nr};
                            Eigen::TensorMap<Eigen::Tensor<Real,1,Eigen::RowMajor>>
                            (fft.getR2C_RMat().data(), quad_nr)=
                            u_spz_BOX_cardinal.slice(start4, count4).reshape(iarray);
                            fft.computeR2C(quad_nr);
                            count4 = {1, quad_nc, 1, 1};
                            iarray = {quad_nc};
                            u_spz_BOX_Fourier.slice(start4, count4).reshape(iarray)=
                            Eigen::TensorMap<Eigen::Tensor<Complex,1,Eigen::RowMajor>>
                            (fft.getR2C_CMat().data(), quad_nc);
                        }
                    }
                }
                
                fft.finalize();
            }
            
            int ncol = 25 * ndim * quad_nc;
            CMatXX out(npts, ncol);
            Eigen::array<int, 2> dims = {npts, ncol};
            Eigen::TensorMap<Eigen::Tensor<Complex,2,Eigen::RowMajor>>
            (out.data(), dims) = u_spz_BOX_Fourier.reshape(dims);
            RMatXX outr = out.real();
            RMatXX outi = out.imag();
            
            
            if (outr.hasNaN() || outi.hasNaN()) {
                throw std::runtime_error("nan in result");
            }
            
            ncw.writeVariableWhole(quad_key + "_RE", outr);
            ncw.writeVariableWhole(quad_key + "_IM", outi);
            
            int iq=std::distance(quad_nrs.begin(), itq);
            std::cout <<  "Done: "<<iq<<" / "<<quad_nrs.size()<<std::endl;
            std::cout.flush();

        }

        ncw.close();
        
        
        exit(0);
        //////////////////////////////////////////
        //////////////////////////////////////////
        //////////////////////////////////////////
        //////////////////////////////////////////
        //////////////////////////////////////////
        //////////////////////////////////////////
        
        
        
        
        
        //////// spectral-element constants
        SpectralConstants::initialize(nPol);  
        
        //////// input parameters 
        int verbose;
        Parameters::buildInparam(pl.mParameters, verbose);
        
        //////// preloop timer
        MultilevelTimer::initialize(Parameters::sOutputDirectory + "/develop/preloop_timer.txt", 4);
        if (pl.mParameters->getValue<bool>("DEVELOP_DIAGNOSE_PRELOOP")) {
            MultilevelTimer::enable();
        }
        
        //////// exodus model and attenuation parameters 
        MultilevelTimer::begin("Build Exodus", 0);
        ExodusModel::buildInparam(pl.mExodusModel, *(pl.mParameters), pl.mAttParameters, verbose);
        MultilevelTimer::end("Build Exodus", 0);
        
        //////// fourier field 
        MultilevelTimer::begin("Build NrField", 0);
        NrField::buildInparam(pl.mNrField, *(pl.mParameters), verbose);
        MultilevelTimer::end("Build NrField", 0);
        
        //////// source
        MultilevelTimer::begin("Build Source", 0);
        Source::buildInparam(pl.mSource, *(pl.mParameters), verbose);
        double srcLat = pl.mSource->getLatitude();
        double srcLon = pl.mSource->getLongitude();
        double srcDep = pl.mSource->getDepth();
        MultilevelTimer::end("Build Source", 0);
        
        //////// 3D models 
        MultilevelTimer::begin("Build 3D Models", 0);
        Volumetric3D::buildInparam(pl.mVolumetric3D, *(pl.mParameters), pl.mExodusModel, 
            srcLat, srcLon, srcDep, verbose);
        Geometric3D::buildInparam(pl.mGeometric3D, *(pl.mParameters), verbose);
        OceanLoad3D::buildInparam(pl.mOceanLoad3D, *(pl.mParameters), verbose);
        MultilevelTimer::end("Build 3D Models", 0);
        
        //////// mesh, phase 1
        // define mesh
        MultilevelTimer::begin("Mesh Definition", 0);
        pl.mMesh = new Mesh(pl.mExodusModel, pl.mNrField, 
            srcLat, srcLon, srcDep, *(pl.mParameters), verbose);
        pl.mMesh->setVolumetric3D(pl.mVolumetric3D);
        pl.mMesh->setGeometric3D(pl.mGeometric3D);
        pl.mMesh->setOceanLoad3D(pl.mOceanLoad3D);
        MultilevelTimer::end("Mesh Definition", 0);
        
        // build unweighted local mesh 
        MultilevelTimer::begin("Build Unweighted Mesh", 0);
        pl.mMesh->buildUnweighted();
        MultilevelTimer::end("Build Unweighted Mesh", 0);
        
        //////// static variables in solver, mainly FFTW
        bool disableWisdomFFTW = pl.mParameters->getValue<bool>("FFTW_DISABLE_WISDOM");
        MultilevelTimer::begin("Initialize FFTW", 0);
        initializeSolverStatic(pl.mMesh->getMaxNr(), disableWisdomFFTW); 
        MultilevelTimer::end("Initialize FFTW", 0);
        
        //////// dt
        MultilevelTimer::begin("Compute DT", 0);
        double dt = pl.mParameters->getValue<double>("TIME_DELTA_T");
        if (dt < tinyDouble) {
            dt = pl.mMesh->getDeltaT();
        }
        double dt_fact = pl.mParameters->getValue<double>("TIME_DELTA_T_FACTOR");
        if (dt_fact < tinyDouble) {
            dt_fact = 1.0;
        }
        dt *= dt_fact;
        MultilevelTimer::end("Compute DT", 0);
        
        XMPI::cout<<dt<<XMPI::endl;
        
        //////// attenuation
        MultilevelTimer::begin("Build Attenuation", 0);
        AttBuilder::buildInparam(pl.mAttBuilder, *(pl.mParameters), pl.mAttParameters, dt, verbose);
        MultilevelTimer::end("Build Attenuation", 0);
        
        //////// mesh, phase 2
        MultilevelTimer::begin("Build Weighted Mesh", 0);
        pl.mMesh->setAttBuilder(pl.mAttBuilder);
        pl.mMesh->buildWeighted();
        MultilevelTimer::end("Build Weighted Mesh", 0);
        
        //////// mesh test 
        // test positive-definiteness and self-adjointness of stiffness and mass matrices
        // better to turn with USE_DOUBLE 
        // pl.mMesh->test();
        // XMPI::barrier();
        // exit(0);
        
        //////// source time function 
        MultilevelTimer::begin("Build Source Time Function", 0);
        STF::buildInparam(pl.mSTF, *(pl.mParameters), dt, verbose);
        MultilevelTimer::end("Build Source Time Function", 0);
        
        //////// receivers
        MultilevelTimer::begin("Build Receivers", 0);
        ReceiverCollection::buildInparam(pl.mReceivers, *(pl.mParameters), 
            srcLat, srcLon, srcDep, pl.mSTF->getSize(), verbose);
        MultilevelTimer::end("Build Receivers", 0);    
        
        //////// computational domain
        MultilevelTimer::begin("Computational Domain", 0);
        sv.mDomain = new Domain();
        
        // release mesh
        MultilevelTimer::begin("Release Mesh", 1);
        pl.mMesh->release(*(sv.mDomain));
        MultilevelTimer::end("Release Mesh", 1);
        
        // release source 
        MultilevelTimer::begin("Release Source", 1);
        pl.mSource->release(*(sv.mDomain), *(pl.mMesh));
        MultilevelTimer::end("Release Source", 1);
        
        // release stf 
        MultilevelTimer::begin("Release STF", 1);
        pl.mSTF->release(*(sv.mDomain));
        MultilevelTimer::end("Release STF", 1);
        
        // release receivers
        MultilevelTimer::begin("Release Receivers", 1);
        pl.mReceivers->release(*(sv.mDomain), *(pl.mMesh), 
            pl.mParameters->getValue<bool>("OUT_STATIONS_DEPTH_REF"));
        MultilevelTimer::begin("Initialize Recorders", 2);
        // sv.mDomain->initializeRecorders();
        MultilevelTimer::end("Initialize Recorders", 2);
        MultilevelTimer::end("Release Receivers", 1);
        
        // verbose domain 
        MultilevelTimer::begin("Verbose", 1);
        if (verbose) {
            XMPI::cout << sv.mDomain->verbose();
        }
        MultilevelTimer::end("Verbose", 1);
        MultilevelTimer::end("Computational Domain", 0);
        
        MultilevelTimer::finalize();
        
        //////////////////////// PREPROCESS DONE ////////////////////////
        
        //////// Newmark
        int infoInt = pl.mParameters->getValue<int>("OPTION_LOOP_INFO_INTERVAL");
        int stabInt = pl.mParameters->getValue<int>("OPTION_STABILITY_INTERVAL");
        bool randomDispl = pl.mParameters->getValue<bool>("DEVELOP_RANDOMIZE_DISP0");
        sv.mNewmark = new Newmark(sv.mDomain, infoInt, stabInt, randomDispl);
        
        ///////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////
        // bin file
        std::stringstream fname;
        #ifdef _SAVE_MEMORY
        if (XMPI::root()) {
            XMPI::mkdir(Parameters::sOutputDirectory + "/bin_save/");
        }
        fname << Parameters::sOutputDirectory << "/bin_save/rank" << XMPI::rank() << ".bin";
        #else 
        if (XMPI::root()) {
            XMPI::mkdir(Parameters::sOutputDirectory + "/bin/");
        }
        fname << Parameters::sOutputDirectory << "/bin/rank" << XMPI::rank() << ".bin";
        #endif
        XMPI::barrier();
        sbs::ofstream ofs(fname.str());
        ofs <<std::string("NewmarkTimeScheme") <<infoInt << stabInt;
        sv.mDomain->dumpToStream(ofs);
        ofs.close();
        
        ///////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////
        
        
        
        
        exit(0);
        //////// final preparations
        // finalize preloop variables before time loop starts
        pl.finalize();
        // forbid matrix allocation in time loop
        #ifndef NDEBUG
            Eigen::internal::set_is_malloc_allowed(false);
        #endif
            
        //////// GoGoGo
        XMPI::barrier();
        sv.mNewmark->solve(verbose);
        
        //////// finalize solver
        // solver 
        sv.mDomain->finalizeRecorders();
        sv.finalize();
        // static variables in solver
        finalizeSolverStatic();
        
        // finalize mpi 
        XMPI::finalize();
        
    } catch (const std::exception &e) {
        // print exception
        XMPI::cout.setp(XMPI::rank());
        XMPI::printException(e);
        
        // abort program
        // TODO 
        // MPI_Abort is necessary here. Otherwise, if an exception
        // is thrown from one of the procs, deadlock will occur.
        // But the problem is, how we free memories on other procs?!
        XMPI::abort();
    }
    
    return 0;
}

#include "SolverFFTW.h"
#include "SolverFFTW_1.h"
#include "SolverFFTW_3.h"
#include "SolverFFTW_N3.h"
#include "SolverFFTW_N6.h"
#include "SolverFFTW_N9.h"
#include "PreloopFFTW.h"
#include "SolidElement.h"
#include "FluidElement.h"

extern void initializeSolverStatic(int maxNr, bool disableWisdomFFTW) {
    // fftw
    SolverFFTW::importWisdom(disableWisdomFFTW);
    SolverFFTW_1::initialize(maxNr);
    SolverFFTW_3::initialize(maxNr); 
    SolverFFTW_N3::initialize(maxNr);
    SolverFFTW_N6::initialize(maxNr);
    SolverFFTW_N9::initialize(maxNr);
    SolverFFTW::exportWisdom();
    // PreloopFFTW::initialize(maxNr);
    // element
    SolidElement::initWorkspace(maxNr / 2);
    FluidElement::initWorkspace(maxNr / 2);
};

extern void finalizeSolverStatic() {
    // fftw
    SolverFFTW_1::finalize();
    SolverFFTW_3::finalize(); 
    SolverFFTW_N3::finalize();
    SolverFFTW_N6::finalize();
    SolverFFTW_N9::finalize();
    PreloopFFTW::finalize();
};

