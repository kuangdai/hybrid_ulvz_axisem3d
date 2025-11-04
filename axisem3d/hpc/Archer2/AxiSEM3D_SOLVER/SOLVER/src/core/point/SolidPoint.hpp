//
//  SolidPoint.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/31/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  solid GLL point

#ifndef SolidPoint_hpp
#define SolidPoint_hpp

#include "Point.hpp"
#include "eigen_element.hpp"

class SolidPoint: public Point {
public:
    // parameter constructor
    SolidPoint(int nr, const eigen::DRow2 &crds, int meshTag,
               std::unique_ptr<const Mass> &mass,
               const TimeScheme &timeScheme);
    
    // stream constructor
    SolidPoint(sbs::ifstream &ifs,
               const TimeScheme &timeScheme);
    
    
    /////////////////////////// mpi ///////////////////////////
    // size for mpi communication
    int sizeComm() const {
        return (int)mFields.mStiff.size();
    }
    
    // feed to mpi buffer
    void feedComm(eigen::CColX &buffer, int &row) const {
        buffer.block(row, 0, mFields.mStiff.size(), 1) =
        Eigen::Map<const eigen::CColX>(mFields.mStiff.data(),
                                       mFields.mStiff.size());
        row += mFields.mStiff.size();
    }
    
    // extract from mpi buffer
    void extractComm(const eigen::CColX &buffer, int &row) {
        mFields.mStiff +=
        Eigen::Map<const eigen::CMatX3>(&buffer(row), mFields.mStiff.rows(), 3);
        row += mFields.mStiff.size();
    }
    
    
    /////////////////////////// time loop ///////////////////////////
    // check stability
    bool stable() const {
        return mFields.mDispl.allFinite();
    }
    
    // stiff to accel
    void computeStiffToAccel();
    
    
    /////////////////////////// element ///////////////////////////
    // scatter displ to element
    void scatterDisplToElement(eigen::vec_ar3_CMatPP_RM &displ,
                               int nu_1_element, int ipol, int jpol) const {
        int nu_1 = mNu + 1;
        
        // copy lower orders
        for (int alpha = 0; alpha < nu_1; alpha++) {
            displ[alpha][0](ipol, jpol) = mFields.mDispl(alpha, 0);
            displ[alpha][1](ipol, jpol) = mFields.mDispl(alpha, 1);
            displ[alpha][2](ipol, jpol) = mFields.mDispl(alpha, 2);
        }
        
        // mask higher orders
        static const numerical::ComplexR czero = 0.;
        for (int alpha = nu_1; alpha < nu_1_element; alpha++) {
            displ[alpha][0](ipol, jpol) = czero;
            displ[alpha][1](ipol, jpol) = czero;
            displ[alpha][2](ipol, jpol) = czero;
        }
    }
    
    // gather stiff from element
    void gatherStiffFromElement(const eigen::vec_ar3_CMatPP_RM &stiff,
                                int ipol, int jpol) {
        // add lower orders only
        for (int alpha = 0; alpha < mNu + 1; alpha++) {
            mFields.mStiff(alpha, 0) -= stiff[alpha][0](ipol, jpol);
            mFields.mStiff(alpha, 1) -= stiff[alpha][1](ipol, jpol);
            mFields.mStiff(alpha, 2) -= stiff[alpha][2](ipol, jpol);
        }
    }
    
    
    /////////////////////////// source ///////////////////////////
    // add force source (external)
    void addForceSource(const eigen::CMatXN3 &force,
                        int nu_1_force, int ipnt) {
        // add minimum orders only
        int nu_1_min = std::min(mNu + 1, nu_1_force);
        mFields.mStiff.topRows(nu_1_min) +=
        force(Eigen::seqN(Eigen::fix<0>, nu_1_min),
              Eigen::seqN(ipnt, Eigen::fix<3>, Eigen::fix<spectral::nPEM>));
    }
    
    // add force source (element)
    void addForceSource(const eigen::vec_ar3_CMatPP_RM &force,
                        int nu_1_force, int ipol, int jpol) {
        // add minimum orders only
        int nu_1_min = std::min(mNu + 1, nu_1_force);
        for (int alpha = 0; alpha < nu_1_min; alpha++) {
            mFields.mStiff(alpha, 0) += force[alpha][0](ipol, jpol);
            mFields.mStiff(alpha, 1) += force[alpha][1](ipol, jpol);
            mFields.mStiff(alpha, 2) += force[alpha][2](ipol, jpol);
        }
    }
    
    
    /////////////////////////// fields ///////////////////////////
    // fields on a solid point
    struct Fields {
        eigen::CMatX3 mStiff = eigen::CMatX3(0, 3);
        eigen::CMatX3 mDispl = eigen::CMatX3(0, 3);
        eigen::CMatX3 mVeloc = eigen::CMatX3(0, 3);
        eigen::CMatX3 mAccel = eigen::CMatX3(0, 3);
    };
    
    // get
    const Fields &getFields() const {
        return mFields;
    }
    
    // set
    Fields &getFields() {
        return mFields;
    }
    
private:
    // fields on a solid point
    Fields mFields;
};

#endif /* SolidPoint_hpp */
