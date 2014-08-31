#/******************************************************************************
#*
#* Blending Liquids : LiquidBlender (LB)
#* Copyright 2014, Karthik Raveendran
#*
#* Cmakefile
#*
#*****************************************************************************/

cmake_minimum_required(VERSION 2.6)
project (liquidBlender)

SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
SET(VERBOSE 1) 

IF(NOT CMAKE_BUILD_TYPE) 
	SET(CMAKE_BUILD_TYPE "Release") 
ELSE()
	MESSAGE("Build-type: '${CMAKE_BUILD_TYPE}'")
ENDIF()


MACRO(ADD_MSVC_PRECOMPILED_HEADER PrecompiledHeader PrecompiledSource SourcesVar)

    GET_FILENAME_COMPONENT(PrecompiledBasename ${PrecompiledHeader} NAME_WE)
    SET(PrecompiledBinary "$(IntDir)/${PrecompiledBasename}.pch")
    SET(Sources ${${SourcesVar}})

    SET_SOURCE_FILES_PROPERTIES(${PrecompiledSource}
                                PROPERTIES COMPILE_FLAGS "/Yc\"${PrecompiledHeader}\" /Fp\"${PrecompiledBinary}\""
                                           OBJECT_OUTPUTS "${PrecompiledBinary}")
    SET_SOURCE_FILES_PROPERTIES(${Sources}
                                PROPERTIES COMPILE_FLAGS "/Yu\"${PrecompiledHeader}\" /FI\"${PrecompiledHeader}\" /Fp\"${PrecompiledBinary}\""
                                           OBJECT_DEPENDS "${PrecompiledBinary}")  
    # Add precompiled header to SourcesVar
    LIST(APPEND ${SourcesVar} ${PrecompiledSource})

ENDMACRO(ADD_MSVC_PRECOMPILED_HEADER)


set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package( Boost REQUIRED COMPONENTS filesystem system)

SET(LB_VC_COMPILE_FLAGS "/MP /W3 /GS- /O2 /Ob2 /Oi /Ot /openmp /fp:precise /D \"_VARIADIC_MAX=10\"")
SET(LB_VC_LINK_FLAGS "") # additional linker flags, not a list

SET(LB_SOURCES		
	src/Animation2D.cpp
	src/Animation3D.cpp
	src/EdgeLabels.cpp
	src/LiquidBlender.cpp
	src/Geometry.cpp
	src/ICP2D.cpp
	src/ICP3D.cpp
	src/MarchingSquares.cpp
	src/Mesh2D.cpp
	src/Mesh3D.cpp
	src/Output.cpp
	#src/stdafx.cpp
	src/Utility.cpp
) 


SET(LB_HEADERS
	src/Animation.h
	src/Blend.hpp
	src/Configuration.hpp
	src/EdgeLabels.h
	src/Geometry.h
	src/Globals.h
	src/ICP.hpp
	src/kdTree.hpp
	src/MarchingSquares2D.h
	src/Mesh.h
	src/Mesh2D.h
	src/Mesh3D.h
	src/nanoflann.hpp
	src/Neighbour.h
	src/Output.hpp
	src/resource.h
	src/Retimer.hpp
	src/Retiming.h
	#src/stdafx.h
	src/targetver.h
	src/Utility.hpp
) 

ADD_MSVC_PRECOMPILED_HEADER("src/stdafx.h" "src/stdafx.cpp" LB_SRCS)

	
# include dirs & libs
SET(LB_INCLUDES
	$ENV{ZLIB_PATH}/include 
	$ENV{EIGEN_PATH}
	$ENV{GOOGLETEST_PATH}/include
	$ENV{GOOGLELOG_PATH}/src/windows	
	${Boost_INCLUDE_DIRS}
	$ENV{LIBJSON_PATH}/include
	$ENV{TBB_PATH}/include
	$ENV{MKL_PATH}/include
	${CMAKE_CURRENT_SOURCE_DIR}/src
)





FIND_LIBRARY(ZLIB zlibwapi $ENV{ZLIB_PATH}/lib/x64)
FIND_LIBRARY(GOOGLELOG libglog $ENV{GOOGLELOG_PATH}/x64/Release)
FIND_LIBRARY(GOOGLETEST gtest $ENV{GOOGLETEST_PATH}/build/Release)
FIND_LIBRARY(GOOGLETEST_MAIN gtest_main $ENV{GOOGLETEST_PATH}/build/Release)
FIND_LIBRARY(JSON libjson $ENV{LIBJSON_PATH}/x64/Release)
FIND_LIBRARY(TBB tbb $ENV{TBB_PATH}/lib/intel64/vc11)
FIND_LIBRARY(MKL_CORE mkl_core_dll $ENV{MKL_PATH}/lib/intel64)
FIND_LIBRARY(MKL_INTEL_THREAD mkl_intel_thread_dll $ENV{MKL_PATH}/lib/intel64)
FIND_LIBRARY(MKL_LP64 mkl_intel_lp64_dll $ENV{MKL_PATH}/lib/intel64)
FIND_LIBRARY(LIBIO libiomp5md.lib $ENV{INTEL_PATH})

SET(LB_LIBS
	${ZLIB}
	${GOOGLELOG}
	${GOOGLETEST}
	${GOOGLETEST_MAIN}
	${JSON}
	${TBB}
	${MKL_CORE}
	${MKL_INTEL_THREAD}
	${MKL_LP64}
	${LIBIO}
	${Boost_FILESYSTEM_LIBRARY}
	${Boost_SYSTEM_LIBRARY}
)


SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${LB_VC_COMPILE_FLAGS}" )
SET(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} ${LB_VC_LINK_FLAGS} ")
SET(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} ${LB_VC_LINK_FLAGS} ")




add_executable( liquidBlender ${LB_SOURCES} )
TARGET_INCLUDE_DIRECTORIES( liquidBlender PUBLIC ${LB_INCLUDES} ) 
TARGET_LINK_LIBRARIES(liquidBlender ${LB_LIBS})



# debug
MESSAGE("DEBUG Flag-Summary - Includes: '${LB_INCLUDES}' | Libs: '${LB_LIBS}' ")
