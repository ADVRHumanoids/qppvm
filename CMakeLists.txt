# CMake minimum version and policies
cmake_minimum_required(VERSION 2.8.12)
if(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 OLD)
    cmake_policy(SET CMP0005 NEW)
    cmake_policy(SET CMP0017 NEW)
endif(COMMAND cmake_policy)

# XBotCore
project(QPPVM_RT_plugin)

# Add cmake dir and include custom macros
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
include(MacroYCMInstallLibrary)

find_package(XCM REQUIRED)
find_package(trajectory_utils REQUIRED)
find_package(OpenSoT REQUIRED)

set(ENABLE_XENO CACHE BOOL "ON")
if ( ${ENABLE_XENO} )
   find_package(Xenomai)
   find_package(RTnet)
endif()

# C++11
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++11" COMPILER_SUPPORTS_CXX11)
check_cxx_compiler_flag("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

# find packages

# include
include_directories(include
                    ${XCM_INCLUDE_DIRS} ${OpenSoT_INCLUDE_DIRS} ${trajectory_utils_INCLUDE_DIRS}
                    )

###########
## Build ##
###########
add_library(QPPVMPlugin SHARED src/QPPVMPlugin.cpp include)
add_library(ForceAccPlugin SHARED src/ForceAcc.cpp include)

if ( Xenomai_FOUND )
#     set_xeno_flags(QPPVMPlugin)
endif()

##########
## Link ##
##########
target_link_libraries(QPPVMPlugin ${XCM_LIBRARIES} ${OpenSoT_LIBRARIES} ${trajectory_utils_LIBRARIES})
target_link_libraries(ForceAccPlugin ${XCM_LIBRARIES} ${OpenSoT_LIBRARIES})



#############
## Install ##
#############

ycm_library_install(QPPVMPlugin 1 0 0)
ycm_library_install(ForceAccPlugin 1 0 0)

include(AddUninstallTarget)
