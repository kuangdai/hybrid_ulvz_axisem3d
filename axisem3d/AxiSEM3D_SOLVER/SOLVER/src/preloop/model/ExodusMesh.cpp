//
//  ExodusMesh.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 9/2/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  Exodus mesh created by salvus mesher

#include "ExodusMesh.hpp"

// tools
#include "timer.hpp"
#include "mpi.hpp"
#include "io.hpp"
#include "bstring.hpp"
#include "eigen_tools.hpp"

// input
#include "InparamYAML.hpp"
#include "NetCDF_Reader.hpp"

// constructor
ExodusMesh::ExodusMesh(const InparamYAML &inparamModel):
mFileName(inparamModel.gets<std::string>("model1D:exodus_mesh")) {
    // memory info
    double memSup = 0.;
    double memAll = 0.;
    
    // opne file on root
    timer::gPreloopTimer.begin("Opening Exodus mesh file");
    NetCDF_Reader reader;
    if (mpi::root()) {
        reader.open(io::gInputDirectory + "/" + mFileName);
    }
    timer::gPreloopTimer.ended("Opening Exodus mesh file");
    
    // read and broadcast
    timer::gPreloopTimer.begin("Reading and broadcasting mesh");
    
    // global variables and records
    timer::gPreloopTimer.begin("Global variables");
    readBcastGlobal(reader, memSup, memAll);
    timer::gPreloopTimer.ended("Global variables");
    
    // connectivity
    timer::gPreloopTimer.begin("Connectivity");
    readBcastConnectivity(reader, memSup, memAll);
    timer::gPreloopTimer.ended("Connectivity");
    
    // coords
    timer::gPreloopTimer.begin("Coordinates");
    readBcastCoordinates(reader, memSup, memAll);
    timer::gPreloopTimer.ended("Coordinates");
    
    // side sets
    timer::gPreloopTimer.begin("Side sets");
    readBcastSideSets(reader, memSup, memAll);
    timer::gPreloopTimer.ended("Side sets");
    
    // elemental variables
    timer::gPreloopTimer.begin("Element variables");
    readBcastElemental(reader, memSup, memAll);
    timer::gPreloopTimer.ended("Element variables");
    
    // radial variables
    timer::gPreloopTimer.begin("Radial variables");
    readBcastRadial(reader, memSup, memAll);
    timer::gPreloopTimer.ended("Radial variables");
    
    // ellipticity
    timer::gPreloopTimer.begin("Ellipticity");
    readBcastEllipticity(reader, memSup, memAll);
    timer::gPreloopTimer.ended("Ellipticity");
    
    // report total memory
    std::stringstream ss;
    ss << "*** Non-scalable memory allocation on super processors: ";
    ss << memSup << " GB ***";
    timer::gPreloopTimer.message(ss.str());
    ss.str("");
    ss << "*** Non-scalable memory allocation on all processors: ";
    ss << memAll << " GB ***";
    timer::gPreloopTimer.message(ss.str());
    
    // finish
    timer::gPreloopTimer.ended("Reading and broadcasting mesh");
}

// verbose
std::string ExodusMesh::verbose() const {
    if (!mpi::root()) {
        // only root has complete data
        return "";
    }
    
    using namespace bstring;
    std::stringstream ss;
    ss << boxTitle("Exodus Mesh");
    
    // overview
    ss << boxSubTitle(0, "Overview");
    ss << boxEquals(2, 20, "Exodus file", mFileName);
    ss << boxEquals(2, 20, "mesh CS", mGlobalRecords.at("crdsys"));
    ss << boxEquals(2, 20, "model name", mGlobalRecords.at("model"));
    ss << boxEquals(2, 20, "isotropic", isIsotropic());
    ss << boxEquals(2, 20, "attenuation", hasAttenuation());
    ss << boxEquals(2, 20, "storage type", mElementNodesStorage ?
                    "element_nodes" : "elements");
    ss << boxEquals(2, 20, "# discontinuities", mDiscontinuities.size());
    
    // geometry
    ss << boxSubTitle(0, "Mesh geometry");
    ss << boxEquals(2, 20, "# elements", getNumQuads());
    ss << boxEquals(2, 20, "# nodes", getNumNodes());
    if (isCartesian()) {
        const std::string &rangeS = range(mNodalCoords.col(0).minCoeff(),
                                          mNodalCoords.col(0).maxCoeff());
        const std::string &rangeZ = range(mNodalCoords.col(1).minCoeff(),
                                          mNodalCoords.col(1).maxCoeff());
        ss << boxEquals(2, 20, "range of s-axis", rangeS);
        ss << boxEquals(2, 20, "range of z-axis", rangeZ);
    } else {
        const auto &r = mNodalCoords.rowwise().norm();
        const auto &t = (mNodalCoords.col(1).cwiseQuotient
                         (r.cwiseMax(numerical::dEpsilon))).array().acos();
        const std::string &rangeR = range(r.minCoeff(), r.maxCoeff());
        const std::string &rangeT = range(t.minCoeff(), t.maxCoeff());
        ss << boxEquals(2, 20, "range of r-axis", rangeR);
        ss << boxEquals(2, 21, "range of θ-axis", rangeT);
    }
    
    // global
    ss << boxSubTitle(0, "Global variables");
    for (auto it = mGlobalVariables.begin();
         it != mGlobalVariables.end(); it++) {
        // hide w and y
        if (it->first.find("w_") != 0 && it->first.find("y_") != 0) {
            ss << boxEquals(2, 20, it->first, it->second);
        }
    }
    for (auto it = mGlobalRecords.begin(); it != mGlobalRecords.end(); it++) {
        ss << boxEquals(2, 20, it->first, it->second);
    }
    
    // element type
    ss << boxSubTitle(0, "Element geometry");
    if (isCartesian()) {
        ss << boxEquals(2, 20, "linear", getNumQuads(), ":");
    } else {
        std::map<std::string, int> eleTypeMap = {
            {"spherical", 0},
            {"linear", 0},
            {"semi-spherical", 0}};
        for (int etype: mElementTypes) {
            if (etype == 0) {
                eleTypeMap["spherical"]++;
            } else if (etype == 1) {
                eleTypeMap["linear"]++;
            } else {
                eleTypeMap["semi-spherical"]++;
            }
        }
        ss << boxEquals(2, 20, eleTypeMap, ":");
    }
    
    // radial
    ss << boxSubTitle(0, "Material properties");
    for (auto it = mRadialVariables.begin();
         it != mRadialVariables.end(); it++) {
        const std::string &prange = range(it->second.minCoeff(),
                                          it->second.maxCoeff());
        ss << boxEquals(2, 20, it->first, prange, "∈");
    }
    
    // side sets
    ss << boxSubTitle(0, "Side sets");
    for (auto it = mSideSets.begin(); it != mSideSets.end(); it++) {
        ss << boxEquals(2, 20, it->first, toString(it->second.size())
                        + " sides", ":");
    }
    
    // miscellaneous (currently only ellipticity)
    ss << boxSubTitle(0, "Miscellaneous");
    ss << boxEquals(2, 20, "ellipticity curve",
                    toString(mEllipticityCurve.cols()) + " knots", ":");
    
    // mesh generation cmdline
    ss << boxBaseline(fmt::gBoxWidth, '-');
    ss << "Mesh generation command line:\n";
    std::vector<std::string> subcmd = split(mCmdMeshGen, " \t");
    std::string line = ">> ";
    std::string tryl;
    for (int icmd = 0; icmd < subcmd.size(); icmd++) {
        tryl = line + subcmd[icmd] + " ";
        if (tryl.size() > fmt::gBoxWidth - 2) {
            ss << line << "\\\n";
            line = "   " + subcmd[icmd] + " ";
        } else {
            line = tryl;
        }
    }
    ss << line << "\n";
    
    // finish
    ss << boxBaseline() << "\n\n";
    return ss.str();
}


///////////////// read and broadcast/////////////////
// global variables
void ExodusMesh::readBcastGlobal(const NetCDF_Reader &reader,
                                 double &memSup, double &memAll) {
    // read
    std::vector<std::string> glbVarNames, glbRecords;
    eigen::DColX glbVarValues;
    if (mpi::root()) {
        reader.readString("name_glo_var", glbVarNames);
        reader.readString("info_records", glbRecords);
        reader.readEigen("vals_glo_var", glbVarValues);
    }
    
    // broadcast
    mpi::bcast(glbVarNames);
    mpi::bcast(glbRecords);
    mpi::bcastEigen(glbVarValues);
    
    // cast variables to map
    for (int igv = 0; igv < glbVarNames.size(); igv++) {
        // "dt" could be misleading
        if (glbVarNames[igv] == "dt") {
            glbVarNames[igv] = "dt (nPol = 1)";
        }
        // make attenuation coefficients readable
        if (glbVarNames[igv].find("w_") == 0 ||
            glbVarNames[igv].find("y_") == 0 ||
            glbVarNames[igv].find("f_") == 0) {
            glbVarNames[igv] = glbVarNames[igv] + "_attenuation";
        }
        mGlobalVariables.insert({glbVarNames[igv], glbVarValues(igv)});
    }
    
    // cast records to map
    const std::vector<std::string> keys = {"crdsys", "model"};
    for (int igr = 0; igr < glbRecords.size(); igr++) {
        // split by "="
        std::vector<std::string> substrs = bstring::split(glbRecords[igr], "=");
        // only include required keys
        if (std::find(keys.begin(), keys.end(), substrs[0]) != keys.end()) {
            mGlobalRecords.insert({substrs[0], substrs[1]});
        }
        // mesh generation cmd
        if (substrs[0].find("cmdl") != std::string::npos) {
            mCmdMeshGen += substrs[1] + " ";
        }
    }
}

// connectivity
void ExodusMesh::readBcastConnectivity(const NetCDF_Reader &reader,
                                       double &memSup, double &memAll) {
    // read
    if (mpi::root()) {
        reader.readEigen("connect1", mConnectivity);
        // let node index start from 0
        mConnectivity.array() -= 1;
    }
    
    // broadcast
    if (mpi::super()) {
        mpi::enterSuper();
        mpi::bcastEigen(mConnectivity);
        mpi::enterWorld();
        // memory info
        timer::gPreloopTimer.message
        (eigen_tools::memoryInfo(mConnectivity, "connectivity", memSup));
    }
}

// coordinates
void ExodusMesh::readBcastCoordinates(const NetCDF_Reader &reader,
                                      double &memSup, double &memAll) {
    // read
    if (mpi::root()) {
        // read
        eigen::DColX x, y;
        reader.readEigen("coordx", x);
        reader.readEigen("coordy", y);
        // cast to matrix
        mNodalCoords.resize(x.rows(), 2);
        mNodalCoords.col(0) = x;
        mNodalCoords.col(1) = y;
    }
    
    // broadcast
    if (mpi::super()) {
        mpi::enterSuper();
        mpi::bcastEigen(mNodalCoords);
        mpi::enterWorld();
        // memory info
        timer::gPreloopTimer.message
        (eigen_tools::memoryInfo(mNodalCoords, "coordinates", memSup));
    }
}

// side sets
void ExodusMesh::readBcastSideSets(const NetCDF_Reader &reader,
                                   double &memSup, double &memAll) {
    // names
    std::vector<std::string> ssNames;
    if (mpi::root()) {
        reader.readString("ss_names", ssNames);
    }
    mpi::bcast(ssNames);
    
    // values
    for (int iss = 0; iss < ssNames.size(); iss++) {
        // read
        eigen::IColX elems, sides;
        if (mpi::root()) {
            const std::string &istr = bstring::toString(iss + 1);
            // elem
            reader.readEigen("elem_ss" + istr, elems);
            // side
            reader.readEigen("side_ss" + istr, sides);
            // let side index start from 0
            elems.array() -= 1;
            sides.array() -= 1;
        }
        
        // broadcast
        mpi::bcastEigen(elems);
        mpi::bcastEigen(sides);
        // memory info
        timer::gPreloopTimer.message
        (eigen_tools::memoryInfo(elems, ssNames[iss], memAll, 2));
        
        // cast to map
        std::map<int, int> mapES;
        for (int ie = 0; ie < elems.rows(); ie++) {
            mapES.insert({elems(ie), sides(ie)});
        }
        mSideSets.insert({ssNames[iss], mapES});
    }
    
    // side set keys
    mKeyLeftSS = isCartesian() ? "x0" : "t0";
    mKeyRightSS = isCartesian() ? "x1" : "t1";
    mKeyBottomSS = isCartesian() ? "y0" : "r0";
    mKeyTopSS = isCartesian() ? "y1" : "r1";
    
    // add minimum edge length to global variable
    double minEdge = std::numeric_limits<double>::max();
    if (mpi::root()) {
        for (int iquad = 0; iquad < getNumQuads(); iquad++) {
            // only consider axial elements
            if (getAxialSide(iquad) >= 0) {
                for (int pA = 0; pA < 4; pA++) {
                    int pB = (pA == 3) ? 0 : pA + 1;
                    const auto &szA =
                    mNodalCoords.row(mConnectivity(iquad, pA));
                    const auto &szB =
                    mNodalCoords.row(mConnectivity(iquad, pB));
                    minEdge = std::min(minEdge, (szA - szB).norm());
                }
            }
        }
    }
    mpi::bcast(minEdge);
    mGlobalVariables.insert({"min_edge_length", minEdge});
}

// elemental variables
void ExodusMesh::readBcastElemental(const NetCDF_Reader &reader,
                                    double &memSup, double &memAll) {
    // only element type
    if (isCartesian()) {
        // all Linear in a Cartesian mesh
        mElementTypes.resize(0);
        return;
    }
    
    // read
    if (mpi::root()) {
        // get exodus-key of "element_type"
        std::vector<std::string> evNames;
        reader.readString("name_elem_var", evNames);
        auto it = std::find(evNames.begin(), evNames.end(), "element_type");
        std::stringstream exKey;
        exKey << "vals_elem_var" << (int)(it - evNames.begin() + 1) << "eb1";
        // read
        eigen::DColX dbuffer;
        reader.readEigen(exKey.str(), dbuffer);
        mElementTypes = dbuffer.cast<int>();
    }
    
    // broadcast
    if (mpi::super()) {
        mpi::enterSuper();
        mpi::bcastEigen(mElementTypes);
        mpi::enterWorld();
        // memory
        timer::gPreloopTimer.message
        (eigen_tools::memoryInfo(mElementTypes, "element_type", memSup));
    }
}

// radial variables
void ExodusMesh::readBcastRadial(const NetCDF_Reader &reader,
                                 double &memSup, double &memAll) {
    //////////////////////// name ////////////////////////
    // names of elemental variables
    std::vector<std::string> evNames;
    if (mpi::root()) {
        reader.readString("name_elem_var", evNames);
    }
    mpi::bcast(evNames);
    
    // name map for quick access
    std::map<std::string, std::string> evNameMap;
    for (int iev = 0; iev < evNames.size(); iev++) {
        std::stringstream exKey;
        exKey << "vals_elem_var" << iev + 1 << "eb1";
        evNameMap.insert({evNames[iev], exKey.str()});
    }
    
    // storage type
    mElementNodesStorage = (evNameMap.find("RHO_0") != evNameMap.end());
    
    // names of radial variables
    std::vector<std::string> rvNames;
    // isIsotropic() is not callable here
    bool isotropic = (evNameMap.find("ETA_0") == evNameMap.end() &&
                      evNameMap.find("ETA") == evNameMap.end());
    if (isotropic) {
        rvNames = std::vector<std::string>({"RHO", "VP", "VS"});
    } else {
        rvNames = std::vector<std::string>({"RHO", "VPV", "VPH",
            "VSV", "VSH", "ETA"});
    }
    if (hasAttenuation()) {
        rvNames.push_back("QMU");
        rvNames.push_back("QKAPPA");
    }
    
    
    //////////////////////// values ////////////////////////
    // dist tolerance of the mesh
    double distTol = mGlobalVariables.at("min_edge_length") / 100.;
    
    // values
    eigen::DMatXX rvValues;
    if (mpi::root()) {
        //////////////////////// depth profile ////////////////////////
        // coords of axial points
        std::vector<double> coords;
        // quad and node tags of axial points
        std::vector<std::pair<int, int>> quad_nodes;
        for (int iquad = 0; iquad < getNumQuads(); iquad++) {
            // only consider axial elements
            int axisSide = getAxialSide(iquad);
            if (axisSide >= 0) {
                // coords
                int otherNode = (axisSide == 3) ? 0 : axisSide + 1;
                double z1 = mNodalCoords(mConnectivity(iquad, axisSide), 1);
                double z2 = mNodalCoords(mConnectivity(iquad, otherNode), 1);
                // negative values not needed
                if (std::min(z1, z2) < 0.) {
                    continue;
                }
                // move both ends inward by tolerance
                z1 += (z2 - z1) / std::abs(z2 - z1) * distTol;
                z2 -= (z2 - z1) / std::abs(z2 - z1) * distTol;
                // add to profile, keep sorted
                auto ir1 = coords.insert
                (std::upper_bound(coords.begin(), coords.end(), z1), z1);
                quad_nodes.insert(quad_nodes.begin() + (ir1 - coords.begin()),
                                  {iquad, axisSide});
                auto ir2 = coords.insert
                (std::upper_bound(coords.begin(), coords.end(), z2), z2);
                quad_nodes.insert(quad_nodes.begin() + (ir2 - coords.begin()),
                                  {iquad, otherNode});
            }
        }
        // cast to matrix
        mRadialCoords = Eigen::Map<eigen::DColX>(coords.data(), coords.size());
        
        
        //////////////////////// depth values ////////////////////////
        rvValues = eigen::DMatXX::Zero(coords.size(), rvNames.size());
        for (int iname = 0; iname < rvNames.size(); iname++) {
            if (mElementNodesStorage) {
                // storage = element_nodes
                std::array<eigen::DColX, 4> buf;
                reader.readEigen(evNameMap.at(rvNames[iname] + "_0"), buf[0]);
                reader.readEigen(evNameMap.at(rvNames[iname] + "_1"), buf[1]);
                reader.readEigen(evNameMap.at(rvNames[iname] + "_2"), buf[2]);
                reader.readEigen(evNameMap.at(rvNames[iname] + "_3"), buf[3]);
                for (int idep = 0; idep < coords.size(); idep++) {
                    rvValues(idep, iname) = buf[quad_nodes[idep].second]
                    (quad_nodes[idep].first);
                }
            } else {
                // storage = elements
                eigen::DColX buf;
                reader.readEigen(evNameMap.at(rvNames[iname]), buf);
                for (int idep = 0; idep < coords.size(); idep++) {
                    rvValues(idep, iname) = buf(quad_nodes[idep].first);
                }
            }
        }
        
        
        //////////////////////// discontinuities ////////////////////////
        // read
        eigen::DColX allDiscs;
        reader.readEigen("discontinuities", allDiscs);
        // to absolute
        if (mGlobalVariables.find("radius") != mGlobalVariables.end()) {
            allDiscs *= mGlobalVariables.at("radius");
        } else {
            allDiscs *= getMeshTop();
        }
        
        // vertical model range
        double min_z, max_z;
        if (isCartesian()) {
            min_z = mNodalCoords.col(1).minCoeff();
            max_z = mNodalCoords.col(1).maxCoeff();
        } else {
            const auto &norm = mNodalCoords.rowwise().norm();
            min_z = norm.minCoeff();
            max_z = norm.maxCoeff();
        }
        
        // remove those out of model range
        std::vector<double> myDiscs;
        for (int idisc = 0; idisc < allDiscs.size(); idisc++) {
            if (allDiscs(idisc) > min_z - distTol &&
                allDiscs(idisc) < max_z + distTol) {
                myDiscs.push_back(allDiscs(idisc));
            }
        }
        mDiscontinuities = Eigen::Map<eigen::DColX>
        (myDiscs.data(), myDiscs.size());
        
        // verify discontinuities
        int disc_found = 0;
        for (int idep = 0; idep < coords.size() - 1; idep++) {
            if (coords[idep + 1] - coords[idep] > distTol * 4) {
                // this gap spans an element
                // do nothing
                continue;
            }
            // check gap between discontinuities
            bool gapIsDisc = false;
            for (int idisc = 0; idisc < mDiscontinuities.size(); idisc++) {
                if (mDiscontinuities[idisc] < coords[idep + 1] &&
                    mDiscontinuities[idisc] > coords[idep]) {
                    // this gap spans a discontinuity
                    // keep it as is
                    gapIsDisc = true;
                    disc_found++;
                    break;
                }
            }
            // this gap is fake
            if (!gapIsDisc) {
                // average the upper and the lower values
                // such difference results from "elements" storage
                rvValues.row(idep) = (rvValues.row(idep) +
                                      rvValues.row(idep + 1)) * .5;
                rvValues.row(idep + 1) = rvValues.row(idep);
            }
        }
        // the top and bottom won't be found
        if (disc_found < mDiscontinuities.size() - 2) {
            throw std::runtime_error("ExodusMesh::readBcastRadial || "
                                     "Error verifying mesh discontinuities.");
        }
    }
    
    // broadcast
    mpi::bcastEigen(mDiscontinuities);
    mpi::bcastEigen(mRadialCoords);
    mpi::bcastEigen(rvValues);
    // memory infor
    timer::gPreloopTimer.message
    (eigen_tools::memoryInfo(mDiscontinuities, "discontinuities", memAll));
    timer::gPreloopTimer.message
    (eigen_tools::memoryInfo(mRadialCoords, "anchoring radii", memAll));
    
    // cast to map
    for (int irv = 0; irv < rvNames.size(); irv++) {
        mRadialVariables.insert({rvNames[irv], rvValues.col(irv)});
        // memory infor
        timer::gPreloopTimer.message
        (eigen_tools::memoryInfo(rvValues.col(irv), rvNames[irv], memAll));
    }
}

// ellipticity
void ExodusMesh::readBcastEllipticity(const NetCDF_Reader &reader,
                                      double &memSup, double &memAll) {
    if (mpi::root()) {
        if (!isCartesian()) {
            // read
            reader.readEigen("ellipticity", mEllipticityCurve);
        } else {
            // constant 1 between [0, 1]
            mEllipticityCurve = eigen::DMatXX_RM::Ones(2, 2);
            mEllipticityCurve(0, 0) = 0.;
        }
        // to absolute
        if (mGlobalVariables.find("radius") != mGlobalVariables.end()) {
            mEllipticityCurve.row(0) *= mGlobalVariables.at("radius");
        } else {
            mEllipticityCurve.row(0) *= getMeshTop();
        }
        // message (broadcast later in geodesy)
        timer::gPreloopTimer.message
        (eigen_tools::memoryInfo(mEllipticityCurve, "ellipticity", memAll));
    }
    mpi::bcastEigen(mEllipticityCurve);
}
