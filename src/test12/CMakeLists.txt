set(CMAKE_CXX_STANDARD 11)

#find_package(EIGEN REQUIRED)
find_package(Eigen3 REQUIRED)
message(STATUS "EIGEN included at ${EIGEN3_INCLUDE_DIR}")


#set(THIRD_PARTY_LIB /home/wasku/Code/cppCode/thirdPartyLibrary)
set(LBFGS_PACKAGE ${THIRD_PARTY_LIB}/LBFGSpp)
set(LBFGS_INCLUDE ${LBFGS_PACKAGE}/include)

include_directories(${LBFGS_INCLUDE})
include_directories(../../include)

#AUX_SOURCE_DIRECTORY(./src/test11 DIR_SRCS)
#AUX_SOURCE_DIRECTORY(./src DIR_SRCS)

add_executable(EigenTest112 main.cpp)
#add_executable(EigenSample src/test/main.cpp ${DIR_SRCS})

