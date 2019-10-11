# Find Spectra library based on the location of GenEigsSolver.h file
#

find_path (SPECTRA_INCLUDE_DIR 
  NAMES Spectra/GenEigsSolver.h
  HINTS "${SPECTRA_ROOT_DIR}/include"
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Spectra DEFAULT_MSG SPECTRA_INCLUDE_DIR)
