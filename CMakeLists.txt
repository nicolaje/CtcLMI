project(Examples)

cmake_minimum_required(VERSION 2.8)


set(CMAKE_CXX_FLAGS "-funroll-all-loops -D_REENTRANT -fPIC -funroll-all-loops")

include_directories(${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/soplex-1.7.2/src ${CMAKE_CURRENT_LIST_DIR}/ibex-dev/include ${CMAKE_CURRENT_LIST_DIR}/sdpa/include ${CMAKE_CURRENT_LIST_DIR}/sdpa/share/sdpa/mumps/build/include ${CMAKE_CURRENT_LIST_DIR}/vibes/)

find_library(IBEX
 NAMES libibex ibex
 PATHS ${CMAKE_CURRENT_LIST_DIR}/ibex-dev/lib
)

find_library(SOPLEX
 NAMES libsoplex soplex
 PATHS ${CMAKE_CURRENT_LIST_DIR}/soplex-1.7.2/lib
)

find_library(CLP
 NAMES libOsiClp Clp libOsiClp.so.1.12.6
 PATHS ${CMAKE_CURRENT_LIST_DIR}/ibex-dev/lib
)

find_library(COINUTILS
 NAMES libCoinUtils CoinUtils
 PATHS ${CMAKE_CURRENT_LIST_DIR}/ibex-dev/lib
)

find_library(SDPA 
 NAMES libsdpa sdpa
 PATHS ${CMAKE_CURRENT_LIST_DIR}/sdpa/lib
)

find_library(MUMPS
 NAMES libdmumps dmumps
 PATHS ${CMAKE_CURRENT_LIST_DIR}/sdpa/share/sdpa/mumps/build/lib/
)

find_library(MUMPS_COMMON
 NAMES libmumps_common mumps_common
 PATHS ${CMAKE_CURRENT_LIST_DIR}/sdpa/share/sdpa/mumps/build/lib/
)

find_library(PORD
 NAMES libpord pord
 PATHS ${CMAKE_CURRENT_LIST_DIR}/sdpa/share/sdpa/mumps/build/lib/
)

find_library(MPISEQ
 NAMES libmpiseq mpiseq
 PATHS ${CMAKE_CURRENT_LIST_DIR}/sdpa/share/sdpa/mumps/build/libseq/
)

find_library(BLAS
 NAMES libopenblas libopenblas_sandybridgep-r0.2.9 openblas
 PATHS ${CMAKE_CURRENT_LIST_DIR}/sdpa/share/sdpa/blas
)

find_library(PRIM
 NAMES libprim prim
 PATHS ${CMAKE_CURRENT_LIST_DIR}/ibex-dev/lib
)

if(SDPA)
 message("SDPA found")
else(SDPA)
 message("SDPA not found")
endif(SDPA)

if(MUMPS)
 message("MUMPS found")
else(MUMPS)
 message("MUMPS not found")
endif(MUMPS)

if(MUMPS_COMMON)
 message("MUMPS_COMMON found")
else(MUMPS_COMMON)
 message("MUMPS_COMMON not found")
endif(MUMPS_COMMON)

if(PORD)
 message("PORD found")
else(PORD)
 message("PORD not found")
endif(PORD)

if(MPISEQ)
 message("MPISEQ found")
else(MPISEQ)
 message("MPISEQ not found")
endif(MPISEQ)

if(BLAS)
 message("BLAS found")
else(BLAS)
 message("BLAS not found")
endif(BLAS)

if(IBEX)
 message("IBEX found")
else(IBEX)
 message("IBEX not found")
endif(IBEX)

if(SOPLEX)
 message("SOPLEX found")
else(SOPLEX)
 message("SOPLEX not found")
endif(SOPLEX)

if(COINUTILS)
 message("COINUTILS found")
else(COINUTILS)
 message("COINUTILS not found")
endif(COINUTILS)

if(CLP)
 message("CLP found")
else(CLP)
 message("CLP not found")
endif(CLP)

if(PRIM)
 message("PRIM found")
else(PRIM)
 message("PRIM not found")
endif(PRIM)

#add_executable(iset iset.cpp ${CMAKE_CURRENT_LIST_DIR}/vibes/vibes.cpp)
add_executable(Ellipses ellipses.cpp CtcBoxLMI.cpp CtcLMI.cpp CtcEllipsoid.cpp CtcEllipsoids.cpp ${CMAKE_CURRENT_LIST_DIR}/vibes/vibes.cpp)

#target_link_libraries(iset ${IBEX} ${SOPLEX} ${SDPA} ${MUMPS} ${MUMPS_COMMON} ${PORD} ${MPISEQ} ${BLAS} ${COINUTILS} ${CLP} ${PRIM} gfortran m quadmath pthread)

target_link_libraries(Ellipses ${IBEX} ${SOPLEX} ${SDPA} ${MUMPS} ${MUMPS_COMMON} ${PORD} ${MPISEQ} ${BLAS} ${COINUTILS} ${CLP} ${PRIM} gfortran m quadmath pthread)

add_executable(Disks disks.cpp CtcLMI.cpp CtcEllipsoid.cpp CtcEllipsoids.cpp ${CMAKE_CURRENT_LIST_DIR}/vibes/vibes.cpp)

target_link_libraries(Disks ${IBEX} ${SOPLEX} ${SDPA} ${MUMPS} ${MUMPS_COMMON} ${PORD} ${MPISEQ} ${BLAS} ${COINUTILS} ${CLP} ${PRIM} gfortran m quadmath pthread)

add_executable(Example1 Example1.cpp CtcBoxLMI.cpp CtcLMI.cpp CtcEllipsoid.cpp CtcEllipsoids.cpp ${CMAKE_CURRENT_LIST_DIR}/vibes/vibes.cpp)

target_link_libraries(Example1 ${IBEX} ${SOPLEX} ${SDPA} ${MUMPS} ${MUMPS_COMMON} ${PORD} ${MPISEQ} ${BLAS} ${COINUTILS} ${CLP} ${PRIM} gfortran m quadmath pthread)

#add_executable(SLAM slam.cpp CtcLMI.cpp CtcEllipsoid.cpp CtcEllipsoids.cpp ${CMAKE_CURRENT_LIST_DIR}/vibes/vibes.cpp)

#target_link_libraries(SLAM ${IBEX} ${SOPLEX} ${SDPA} ${MUMPS} ${MUMPS_COMMON} ${PORD} ${MPISEQ} ${BLAS} ${COINUTILS} ${CLP} ${PRIM} gfortran m quadmath pthread)

