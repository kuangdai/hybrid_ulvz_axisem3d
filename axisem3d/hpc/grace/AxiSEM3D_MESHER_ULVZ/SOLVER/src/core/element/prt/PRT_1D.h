// PRT_1D.h
// created by Kuangdai on 19-May-2017 
// 1D particle relabelling transformation 

#pragma once

#include "PRT.h"
#include "eigenc.h"

class PRT_1D: public PRT {
public:
    PRT_1D(const std::array<RMatPP, 5> &X): mXStruct(X) {};
    ~PRT_1D() {};
    
    void sphericalToUndulated(FluidResponse &response) const;
    void undulatedToSpherical(FluidResponse &response) const;
    
    void sphericalToUndulated(SolidResponse &response) const;
    void undulatedToSpherical(SolidResponse &response) const;
    
    // 9 to 9, for curl computation
    void sphericalToUndulated9(SolidResponse &response) const;
    
    std::string verbose() const {return "PRT_1D";};
    bool is1D() const {return true;};
    
    void dumpToStream(sbs::ofstream &ofs) const {
        RMatXN zero(0, nPntElem);
        ofs<<true;
        
        ofs<<true;
        ofs<<mXStruct[0];
        ofs<<zero;
        
        ofs<<true;
        ofs<<mXStruct[1];
        ofs<<zero;
        
        ofs<<true;
        ofs<<mXStruct[2];
        ofs<<zero;
        
        ofs<<true;
        ofs<<mXStruct[3];
        ofs<<zero;
        
        #if _SAVE_MEMORY
        
        ofs<<true;
        ofs<<mXStruct[4];
        ofs<<zero;
        
        #else
        
        ofs<<true;
        ofs<<(mXStruct[0].cwiseProduct(mXStruct[4])).eval();
        ofs<<zero;
        
        ofs<<true;
        ofs<<(mXStruct[1].cwiseProduct(mXStruct[4])).eval();
        ofs<<zero;
        
        ofs<<true;
        ofs<<(mXStruct[2].cwiseProduct(mXStruct[4])).eval();
        ofs<<zero;
        
        ofs<<true;
        ofs<<(mXStruct[3].cwiseProduct(mXStruct[4])).eval();
        ofs<<zero;
        
        #endif
    }

private:
    std::array<RMatPP, 5> mXStruct;
};

