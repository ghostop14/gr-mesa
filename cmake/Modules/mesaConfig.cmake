INCLUDE(FindPkgConfig)
PKG_CHECK_MODULES(PC_MESA mesa)

FIND_PATH(
    MESA_INCLUDE_DIRS
    NAMES mesa/api.h
    HINTS $ENV{MESA_DIR}/include
        ${PC_MESA_INCLUDEDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/include
          /usr/local/include
          /usr/include
)

FIND_LIBRARY(
    MESA_LIBRARIES
    NAMES gnuradio-mesa
    HINTS $ENV{MESA_DIR}/lib
        ${PC_MESA_LIBDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/lib
          ${CMAKE_INSTALL_PREFIX}/lib64
          /usr/local/lib
          /usr/local/lib64
          /usr/lib
          /usr/lib64
          )

include("${CMAKE_CURRENT_LIST_DIR}/mesaTarget.cmake")

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(MESA DEFAULT_MSG MESA_LIBRARIES MESA_INCLUDE_DIRS)
MARK_AS_ADVANCED(MESA_LIBRARIES MESA_INCLUDE_DIRS)
