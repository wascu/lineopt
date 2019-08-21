cmake_minimum_required(VERSION 3.1)
project(EigenSample)

set(CMAKE_CXX_STANDARD 11)



#
#find_package(EIGEN REQUIRED)
#message(STATUS "EIGEN included at ${EIGEN3_INCLUDE_DIR}")


set(THIRD_PARTY_LIB /home/wasku/Code/cppCode/thirdPartyLibrary)
set(LBFGS_PACKAGE ${THIRD_PARTY_LIB}/LBFGSpp)
set(LBFGS_INCLUDE ${LBFGS_PACKAGE}/include)

include_directories(${LBFGS_INCLUDE})

add_executable(EigenSample main.cpp src/TCallback.cpp include/TCallback.h)

