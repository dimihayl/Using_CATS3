CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
#It is recommended that you copy all CMake related files in a local folder of yours,
#else whenever you pull/push there will be conflicts to resolve

SET(PROJECT_NAME "dlm_hst_cnv")

project(${PROJECT_NAME})
# SET PATHS #
SET(PROJECT_DESTINATION "${CMAKE_CURRENT_LIST_DIR}/../")
SET(PROJECT_USINGCATS3 "${CMAKE_CURRENT_LIST_DIR}/../")#the destination of your project
SET(CATS_PATH $ENV{CATS3})
#SET(CATS_PATH "/home/dimihayl/Software/Repositories/DLM/install")

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_DESTINATION}/bin)

add_executable(${PROJECT_NAME} ${PROJECT_DESTINATION}/dlm_hst_cnv.cpp
#ADD BELOW ANY OTHER .CPP FILES THAT YOU WOULD LIKE TO COMPILE (FULL PATH)
)
# DO NOT CHANGE THE REST #

file(STRINGS ${CATS_PATH}/bin/CMakeDLM.txt CMakeDLM)
list(GET CMakeDLM 2 ROOT_PATH)
message(STATUS "Using this ROOT location: ${ROOT_PATH}")

list(GET CMakeDLM 3 GSL_INCLUDE)
message(STATUS "Using this GSL include location: ${GSL_INCLUDE}")

list(GET CMakeDLM 4 GSL_LIB)
message(STATUS "Using this GSL lib location: ${GSL_LIB}")

execute_process(COMMAND bash -c "${ROOT_PATH}/bin/root-config --cflags" OUTPUT_VARIABLE CFLAGS)
execute_process(COMMAND bash -c "${ROOT_PATH}/bin/root-config --libs" OUTPUT_VARIABLE LIBS)
execute_process(COMMAND bash -c "${ROOT_PATH}/bin/root-config --glibs" OUTPUT_VARIABLE GLIBS)
execute_process(COMMAND bash -c "${ROOT_PATH}/bin/root-config --incdir" OUTPUT_VARIABLE ROOT_INCLUDE)
execute_process(COMMAND bash -c "${CATS_PATH}/bin/cats-config --incdir" OUTPUT_VARIABLE CATS_INCLUDE)
execute_process(COMMAND bash -c "${CATS_PATH}/bin/cats-config --libs" OUTPUT_VARIABLE CATS_LIBS)

string(REGEX REPLACE "\n$" "" CFLAGS "${CFLAGS}")
string(REGEX REPLACE "\n$" "" LIBS "${LIBS}")
string(REGEX REPLACE "\n$" "" GLIBS "${GLIBS}")
string(REGEX REPLACE "\n$" "" ROOT_INCLUDE "${ROOT_INCLUDE}")
string(REGEX REPLACE "\n$" "" CATS_INCLUDE "${CATS_INCLUDE}")
string(REGEX REPLACE "\n$" "" CATS_LIBS "${CATS_LIBS}")

#string(APPEND CFLAGS " -O2 -std=c++11 -lgomp -pthread -fopenmp  ")
set(CFLAGS " -O2 -std=c++11 -pthread ${CFLAGS}")
#string(APPEND CFLAGS " -O2 -std=c++11 -pthread ")

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS ${CFLAGS})

SET(VERSION_MAJOR 1)
SET(VERSION_MINOR 0)
SET(VERSION_PATCH 0)
SET(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR})

include_directories(${PROJECT_USINGCATS3})
include_directories(${PROJECT_DESTINATION})
include_directories(${CATS_INCLUDE})
include_directories(${GSL_INCLUDE})
include_directories(${ROOT_INCLUDE})

target_link_libraries(${PROJECT_NAME} -L${CATS_LIBS} ${LIBS} -L${GSL_LIB} -lgsl -lgslcblas)
