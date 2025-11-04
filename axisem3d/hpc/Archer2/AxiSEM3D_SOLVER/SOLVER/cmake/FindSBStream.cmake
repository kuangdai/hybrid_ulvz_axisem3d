# Kuangdai created this

# - Find the SimpleBinStream library
#
# Usage:
#   find_package(SBStream [REQUIRED] [QUIET])
#   using SBSTREAM_ROOT to hint the path
#
# It sets the following variables:
#   SBStream_INCLUDE_DIR  ... SimpleBinStream include directory


# only need include dir
if(NOT SBStream_INCLUDE_DIR)
find_path(
SBStream_INCLUDE_DIR
NAMES SimpleBinStream.h
HINTS
${SBSTREAM_ROOT}
$ENV{SBSTREAM_ROOT}
PATH_SUFFIXES TestBinStream
)
endif()

# standard handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SimpleBinStream DEFAULT_MSG
SBStream_INCLUDE_DIR)
mark_as_advanced(SBStream_INCLUDE_DIR)
