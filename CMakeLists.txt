cmake_minimum_required(VERSION 3.17)
project(digitizationpp VERSION 0.0.1)

find_package(ROOT REQUIRED)
find_package(OpenCV CONFIG REQUIRED)
find_package(TBB CONFIG REQUIRED)

add_compile_options(-Wall -Wextra -O3)

#Finds c++ standard from ROOT and forces it to the compiler
if(DEFINED ROOT_CXX_FLAGS)
    #Extract  c++std info
    string(REGEX MATCH "-std=c\\+\\+[0-9]+" STD_FLAGS ${ROOT_CXX_FLAGS})
    string(REGEX REPLACE "-std=c\\+\\+" "" ROOT_CXX ${STD_FLAGS})
endif()

if(ROOT_CXX STREQUAL "20")
    set(CMAKE_CXX_STANDARD 20)
elseif(ROOT_CXX STREQUAL "17")
    set(CMAKE_CXX_STANDARD 17)
else()
    message(FATAL_ERROR "ROOT compiled with non-supported C++ version. c++17 and 20 are supported")
endif()

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -no-pie")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -no-pie") 

include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR} ${OpenCV_INCLUDE_DIRS} $ENV{ROOTANASYS}/include)

#Debug comments
#message(ROOT_CXX=${ROOT_CXX})
#message(ROOTANASYS="$ENV{ROOTANASYS}/include")
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    if (_variableName MATCHES "ROOT")
#        message(STATUS "${_variableName}=${${_variableName}}")
#    endif()
#endforeach()

link_directories(${ROOT_LIBRARY_DIR}
                 "$ENV{ROOTANASYS}/lib"
                )

# --------    s3  --------
add_library(s3 "${PROJECT_SOURCE_DIR}/src/s3.cxx")
target_include_directories(s3 PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          "${PROJECT_SOURCE_DIR}/include"
                          )
                                                    
target_link_libraries(s3 PUBLIC curl)
                          
# -------- cygnolib --------
add_library(cygnolib
           ${PROJECT_SOURCE_DIR}/src/cygnolib.cxx
           )
target_include_directories(cygnolib PUBLIC
                          ${PROJECT_BINARY_DIR}
                          ${PROJECT_SOURCE_DIR}/include
                          $ENV{ROOTANASYS}/include
                          )
target_link_libraries(cygnolib PUBLIC ${ROOT_LIBRARIES} rootana z s3)


add_executable(digitizationpp
               "${PROJECT_SOURCE_DIR}/main.cxx"
               "${PROJECT_SOURCE_DIR}/src/ConfigManager.cxx"
               "${PROJECT_SOURCE_DIR}/src/DigitizationRunner.cxx"
               "${PROJECT_SOURCE_DIR}/src/Globals.cxx"#
               "${PROJECT_SOURCE_DIR}/src/TrackProcessor.cxx"
               "${PROJECT_SOURCE_DIR}/src/Utils.cxx"
               )
target_link_libraries(digitizationpp PUBLIC ${ROOT_LIBRARIES} ${OpenCV_LIBS} s3 curl rootana cygnolib z opencv_imgcodecs opencv_core tbb)

add_executable(joinPeds "${PROJECT_SOURCE_DIR}/joinPeds.cxx")
target_link_libraries(joinPeds PUBLIC ${ROOT_LIBRARIES} ${OpenCV_LIBS} s3 curl rootana cygnolib z opencv_imgcodecs opencv_core)

