//
//  io.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/8/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  IO interfaces

#ifndef io_hpp
#define io_hpp

#include <string>
#include <stdexcept>

namespace io {
    //////////////////// input/output dirs ////////////////////
    extern std::string gInputDirectory;
    extern std::string gOutputDirectory;
    
    // check dir existence
    bool dirExists(const std::string &path);
    
    // mkdir
    void mkdir(const std::string &path);
    
    // verify input/output dirs under executable dir
    void verifyDirectories(int argc, char *argv[], std::string &warning);
    
    
    //////////////////// source code dir ////////////////////
    const std::string gProjectDirectory = _PROJ_DIR;
    
    
    //////////////////// runtime verbose ////////////////////
    // verbose levels
    enum class VerboseLevel {
        none = 0,
        essential = 1,
        detailed = 2
    };
    
    // verbose control
    extern VerboseLevel gVerbose;
    extern bool gVerboseWarnings;
    
    // welcome message
    std::string welcome();
    
    // verbose
    std::string verbose();
    
    
    //////////////////// cout on root ////////////////////
    // this allows root-only printing with
    // io::cout << "something" << io::endl;
    struct mpi_root_cout {
        // constructor
        mpi_root_cout();
        
        // mimic cout
        template <typename Type>
        const mpi_root_cout &operator<<(const Type &val) const {
            if (mMyWorldRank == mCoutWorldRank) {
                (*mCoutStream) << val;
            }
            return *this;
        }
        
        // change cout rank
        void setMyWorldRank(int rank) {
            if (mMyWorldRank < 0) {
                mMyWorldRank = rank;
            } else {
                throw std::runtime_error("mpi_root_cout::setMyWorldRank || "
                                         "World rank can be set only once.");
            }
        }
        
        // change cout rank
        void setCoutWorldRank(int rank) {
            mCoutWorldRank = rank;
        }
        
    private:
        // my world rank
        int mMyWorldRank = -1;
        
        // world rank for cout
        int mCoutWorldRank = 0;
        
        // stream: avoid referring std::cout in .hpp
        // because <iostream> is a heavy header
        std::ostream *mCoutStream = 0;
    };
    
    extern mpi_root_cout cout;
    const std::string endl = "\n";
}

#endif /* io_hpp */
