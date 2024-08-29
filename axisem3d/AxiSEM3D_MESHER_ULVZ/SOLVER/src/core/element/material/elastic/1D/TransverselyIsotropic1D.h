// TransverselyIsotropic1D.h
// created by Kuangdai on 22-Apr-2016 
// transversely isotropic 1D material

#pragma once

#include "Elastic1D.h"
#include "eigenc.h"

class TransverselyIsotropic1D: public Elastic1D {
public:
    // constructor
    TransverselyIsotropic1D(const RMatPP &A, const RMatPP &C, const RMatPP &F, 
        const RMatPP &L, const RMatPP &N, Attenuation1D *att):
        Elastic1D(att) , mA(A), mC(C), mF(F), mL(L), mN(N), mN2(two * N) {};
    
    // STEP 2: strain ==>>> stress
    void strainToStress(SolidResponse &response) const;
    
    // verbose
    std::string verbose() const {return "TransverselyIsotropic1D";};
    
    // need TIso
    bool needTIso() const {return true;};
    
    void dumpToStream(sbs::ofstream &ofs,const RRowN &ifact) const  {
        ofs << "TransverselyIsotropic";
        
        Elastic1D::dumpToStream(ofs, ifact);
            
        RMatPP temp;
        int ipnt;
        ipnt = 0;
        for (int ipol = 0; ipol <= nPol; ipol++) {
            for (int jpol = 0; jpol <= nPol; jpol++) {
                temp(ipol,jpol) = mA(ipol, jpol);
                ipnt++;
            }
        }
        ofs << true<<temp<<0<<25;
        
        ipnt = 0;
        for (int ipol = 0; ipol <= nPol; ipol++) {
            for (int jpol = 0; jpol <= nPol; jpol++) {
                temp(ipol,jpol) = mC(ipol, jpol);
                ipnt++;
            }
        }
        ofs << true<<temp<<0<<25;
        
        ipnt = 0;
        for (int ipol = 0; ipol <= nPol; ipol++) {
            for (int jpol = 0; jpol <= nPol; jpol++) {
                temp(ipol,jpol) = mF(ipol, jpol) ;
                ipnt++;
            }
        }
        ofs << true<<temp<<0<<25;
        
        ipnt = 0;
        for (int ipol = 0; ipol <= nPol; ipol++) {
            for (int jpol = 0; jpol <= nPol; jpol++) {
                temp(ipol,jpol) = mL(ipol, jpol);
                ipnt++;
            }
        }
        ofs << true<<temp<<0<<25;
        
        ipnt = 0;
        for (int ipol = 0; ipol <= nPol; ipol++) {
            for (int jpol = 0; jpol <= nPol; jpol++) {
                temp(ipol,jpol) = mN(ipol, jpol);
                ipnt++;
            }
        }
        ofs << true<<temp<<0<<25;
        
        #ifndef _SAVE_MEMORY
        ipnt = 0;
        for (int ipol = 0; ipol <= nPol; ipol++) {
            for (int jpol = 0; jpol <= nPol; jpol++) {
                temp(ipol,jpol) =(- mN2(ipol, jpol));
                ipnt++;
            }
        }
        ofs << true<<temp<<0<<25;
        #endif
        
        
    }
    
private:
    
    // Cijkl scaled by integral factor
    RMatPP mA;
    RMatPP mC;
    RMatPP mF;
    RMatPP mL;
    RMatPP mN;
    RMatPP mN2;
};
