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
add_library(OpenSotIkTestPlugin SHARED src/OpenSotIkTestPlugin.cpp include)

if ( Xenomai_FOUND )
#     set_xeno_flags(QPPVMPlugin)
endif()

##########
## Link ##
##########
target_link_libraries(QPPVMPlugin ${XCM_LIBRARIES} ${OpenSoT_LIBRARIES} ${trajectory_utils_LIBRARIES})
target_link_libraries(OpenSotIkTestPlugin ${XCM_LIBRARIES} ${OpenSoT_LIBRARIES})

#############
## Install ##
#############

# QPPVM PLUGIN ##################################################

set(VARS_PREFIX "QPPVMPlugin")

set(QPPVMPlugin_MAJOR_VERSION 0)
set(QPPVMPlugin_MINOR_VERSION 1)
set(QPPVMPlugin_PATCH_VERSION 0)
set(QPPVMPlugin_VERSION ${QPPVMPlugin_MAJOR_VERSION}.${QPPVMPlugin_MINOR_VERSION}.${QPPVMPlugin_PATCH_VERSION})

find_package(YCM REQUIRED)
include(YCMDefaultDirs)
ycm_default_dirs(${VARS_PREFIX})

target_include_directories(QPPVMPlugin PUBLIC "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>"
                                                "$<INSTALL_INTERFACE:${${VARS_PREFIX}_INSTALL_INCLUDEDIR}>")

set_target_properties(QPPVMPlugin PROPERTIES VERSION ${${VARS_PREFIX}_VERSION}
                                               SOVERSION ${${VARS_PREFIX}_VERSION})

install(DIRECTORY include/
        DESTINATION "${${VARS_PREFIX}_INSTALL_INCLUDEDIR}"
        FILES_MATCHING PATTERN "*.h*")

install(TARGETS QPPVMPlugin
        EXPORT QPPVMPlugin
        ARCHIVE DESTINATION "${${VARS_PREFIX}_INSTALL_BINDIR}" COMPONENT lib
        RUNTIME DESTINATION "${${VARS_PREFIX}_INSTALL_BINDIR}" COMPONENT bin
        LIBRARY DESTINATION "${${VARS_PREFIX}_INSTALL_LIBDIR}" COMPONENT shlib)



#enabling it will add all QPPVMPlugin dependencies as dependencies for third party users
set_property(GLOBAL APPEND PROPERTY ${VARS_PREFIX}_TARGETS QPPVMPlugin)


include(InstallBasicPackageFiles)
install_basic_package_files(QPPVMPlugin VARS_PREFIX ${VARS_PREFIX}
                                    VERSION ${${VARS_PREFIX}_VERSION}
                                    COMPATIBILITY SameMajorVersion
                                    TARGETS_PROPERTY ${VARS_PREFIX}_TARGETS
                                    NO_CHECK_REQUIRED_COMPONENTS_MACRO)

# OPENSOTIK PLUGIN ##################################################


set(VARS_PREFIX "OpenSotIkTestPlugin")

set(OpenSotIkTestPlugin_MAJOR_VERSION 0)
set(OpenSotIkTestPlugin_MINOR_VERSION 1)
set(OpenSotIkTestPlugin_PATCH_VERSION 0)
set(OpenSotIkTestPlugin_VERSION ${OpenSotIkTestPlugin_MAJOR_VERSION}.${OpenSotIkTestPlugin_MINOR_VERSION}.${OpenSotIkTestPlugin_PATCH_VERSION})

find_package(YCM REQUIRED)
include(YCMDefaultDirs)
ycm_default_dirs(${VARS_PREFIX})

target_include_directories(OpenSotIkTestPlugin PUBLIC "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>"
                                                "$<INSTALL_INTERFACE:${${VARS_PREFIX}_INSTALL_INCLUDEDIR}>")

set_target_properties(OpenSotIkTestPlugin PROPERTIES VERSION ${${VARS_PREFIX}_VERSION}
                                               SOVERSION ${${VARS_PREFIX}_VERSION})

install(DIRECTORY include/
        DESTINATION "${${VARS_PREFIX}_INSTALL_INCLUDEDIR}"
        FILES_MATCHING PATTERN "*.h*")

install(TARGETS OpenSotIkTestPlugin
        EXPORT OpenSotIkTestPlugin
        ARCHIVE DESTINATION "${${VARS_PREFIX}_INSTALL_BINDIR}" COMPONENT lib
        RUNTIME DESTINATION "${${VARS_PREFIX}_INSTALL_BINDIR}" COMPONENT bin
        LIBRARY DESTINATION "${${VARS_PREFIX}_INSTALL_LIBDIR}" COMPONENT shlib)



#enabling it will add all OpenSotIkTestPlugin dependencies as dependencies for third party users
set_property(GLOBAL APPEND PROPERTY ${VARS_PREFIX}_TARGETS OpenSotIkTestPlugin)


include(InstallBasicPackageFiles)
install_basic_package_files(OpenSotIkTestPlugin VARS_PREFIX ${VARS_PREFIX}
                                    VERSION ${${VARS_PREFIX}_VERSION}
                                    COMPATIBILITY SameMajorVersion
                                    TARGETS_PROPERTY ${VARS_PREFIX}_TARGETS
                                    NO_CHECK_REQUIRED_COMPONENTS_MACRO)
include(AddUninstallTarget)
