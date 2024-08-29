# - Find the SimpleBinStream library
#
# Usage:
#   find_package(SimpleBinStream [REQUIRED] [QUIET] )
#     
# It sets the following variables:
#   SimpleBinStream_FOUND               ... true if SimpleBinStream is found on the system
#   SimpleBinStream_INCLUDE_DIR         ... SimpleBinStream include directory

if(NOT SimpleBinStream_INCLUDE_DIR)
    
    #find includes
    find_path(
        SimpleBinStream_INCLUDE_DIR
        NAMES SimpleBinStream.h
        HINTS 
        ${SimpleBinStream_ROOT}
        $ENV{SimpleBinStream_ROOT}
        PATH_SUFFIXES TestBinStream
    )
    
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SimpleBinStream DEFAULT_MSG
SimpleBinStream_INCLUDE_DIR)

mark_as_advanced(SimpleBinStream_INCLUDE_DIR)

