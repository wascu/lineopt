cmake_minimum_required(VERSION 3.1)
#project(EigenSample)
set(CMAKE_CXX_STANDARD 11)

find_package(Eigen3 REQUIRED)
message(STATUS "EIGEN included at ${EIGEN3_INCLUDE_DIR}")

#set(THIRD_PARTY_LIB /home/wasku/Code/cppCode/thirdPartyLibrary)
set(LBFGS_PACKAGE ${THIRD_PARTY_LIB}/LBFGSpp)
set(LBFGS_INCLUDE ${LBFGS_PACKAGE}/include)


message(STATUS "LBFGS included at ${LBFGS_INCLUDE}")
include_directories(${LBFGS_INCLUDE})
include_directories(../../include)

add_executable(eigenTest11 main.cpp Jfunc.cpp MdisFunc.cpp)
target_link_libraries(eigenTest11 Eigen3::Eigen)