PROJECT(stap++)

cmake_minimum_required(VERSION 3.0)

OPTION(STAP++_DEBUG "Output detailed results for debugging." OFF)
IF(STAP++_DEBUG)
   ADD_DEFINITIONS(-D_DEBUG_)
ENDIF()

OPTION(STAP++_TEST "Output test infomation to validate." OFF)
IF(STAP++_TEST)
    ADD_DEFINITIONS(-D__TEST__)
ENDIF()

OPTION(USE_MKL "Use MKL PARDISO solver instead of LDLT solver." OFF)
IF(USE_MKL)
    ADD_DEFINITIONS(-DMKL)
    IF(UNIX)
        INCLUDE_DIRECTORIES(/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/include)
        LINK_DIRECTORIES(/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64)
    ELSE()
        INCLUDE_DIRECTORIES("C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2018.1.156/windows/mkl/include")
        LINK_DIRECTORIES("C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2018.1.156/windows/mkl/lib/intel64_win")
    ENDIF()
ENDIF()

ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
ADD_COMPILE_OPTIONS(-std=c++11)

INCLUDE_DIRECTORIES(eigen-3)
INCLUDE_DIRECTORIES(h)

AUX_SOURCE_DIRECTORY(cpp SRC)
AUX_SOURCE_DIRECTORY(cpp/Elements SRC_ELEMENTS)
FILE(GLOB_RECURSE HEAD h/*.h)

source_group(SOURCE\ FILES FILES ${SRC})
source_group(ELEMENTS\ SOURCE\ FILES FILES ${SRC_ELEMENTS})
source_group(HEADER\ FILES FILES ${HEAD})

ADD_EXECUTABLE(stap++ ${SRC} ${SRC_ELEMENTS} ${HEAD})

IF(USE_MKL)
    TARGET_LINK_LIBRARIES(stap++ mkl_intel_lp64 mkl_sequential mkl_core)
ENDIF()
