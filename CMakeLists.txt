cmake_minimum_required(VERSION 2.6)
project(main)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
#find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net RooStats)
if(ROOT_FOUND)
        message("ROOT FOUND")
        include(${ROOT_USE_FILE})
        include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
        add_definitions(${ROOT_CXX_FLAGS})
endif()


INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)

add_library(mylib SHARED ${sources})
target_link_libraries(mylib ${ROOT_LIBRARIES})

add_executable(main main.C ${sources})
target_link_libraries(main mylib)
