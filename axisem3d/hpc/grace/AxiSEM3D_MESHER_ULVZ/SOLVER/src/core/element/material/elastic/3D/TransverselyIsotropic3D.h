// TransverselyIsotropic3D.h
// created by Kuangdai on 22-Apr-2016 
// transversely isotropic 3D material 

#pragma once

#include "Elastic3D.h"
#include "eigenc.h"

class TransverselyIsotropic3D: public Elastic3D {
public:
    // constructor
    TransverselyIsotropic3D(const RMatXN &A, const RMatXN &C, const RMatXN &F, 
        const RMatXN &L, const RMatXN &N, Attenuation3D *att):
        Elastic3D(att) , mA(A), mC(C), mF(F), mL(L), mN(N), mN2(two * N) {};
        
    // STEP 2: strain ==>>> stress
    void strainToStress(SolidResponse &response) const;
    
    // check compatibility
    void checkCompatibility(int Nr) const; 
    
    // verbose
    std::string verbose() const {return "TransverselyIsotropic3D";};
    
    // need TIso
    bool needTIso() const {return true;};
    
    void dumpToStream(sbs::ofstream &ofs,const RRowN &ifact) const  {
        ofs << "TransverselyIsotropic";
        Elastic3D::dumpToStream(ofs, ifact);
        RMatXN a = mA;
        RMatXN c = mC;
        RMatXN f = mF;
        RMatXN l = mL;
        RMatXN n = mN;
        RMatXN n2 = -mN2;
        
        // a.array().rowwise()/=ifact.array();
        // c.array().rowwise()/=ifact.array();
        // f.array().rowwise()/=ifact.array();
        // l.array().rowwise()/=ifact.array();
        // n.array().rowwise()/=ifact.array();
        // n2.array().rowwise()/=ifact.array();
        // 
        ofs<<false<<a;
        ofs<<false<<c;
        ofs<<false<<f;
        ofs<<false<<l;
        ofs<<false<<n;
        
        #ifndef _SAVE_MEMORY
        ofs<<false<<n2;
        #endif
    }
    
private:
    // Cijkl scaled by integral factor
    RMatXN mA;
    RMatXN mC;
    RMatXN mF;
    RMatXN mL;
    RMatXN mN;
    RMatXN mN2;
};
