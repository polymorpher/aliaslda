cmake_minimum_required(VERSION 3.1)
project(aliaslda)
include_directories(/usr/local/include)
include_directories(/usr/local/lib)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#set(CMAKE_BUILD_TYPE Release)
set(CMAKE_BUILD_TYPE)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

set(SOURCE_FILES
    #AliasHDPLDA.hpp
    #AliasLDA_Space.hpp
    #AliasLDA_Speed.hpp
    AliasLDA_UltraSpeed.hpp
    AliasLDATester.cpp
    AliasTable.cpp
    AliasTable.hpp
    #AliasWordPDPLDA.cpp
    #AliasWordPDPLDA.hpp
    DataGroup.hpp
    Dictionary.hpp
    FileReader.cpp
    FileReader.hpp
    #FlatMapTest.cpp
    #GeneralizedStirlingNumber.hpp
    #HDPLDA.hpp
    ILDA.hpp
    #LDA.cpp
    LDA.hpp
    SparseLDA_Speed.hpp
    SparseLDA_Speed.cpp
    Logger.hpp
    PerplexityAnalyzer.hpp
    #SparseLDA_Space.cpp
    #SparseLDA_Space.hpp
    #SparseLDA_Speed.cpp
    #SparseLDA_Speed.hpp
    #StirlingNumberTester.cpp
    #StructTester.hpp
    util.hpp
    #WordPDPLDA.cpp
    #WordPDPLDA.hpp
    )
find_package(JsonCpp REQUIRED)
find_package(Boost REQUIRED regex)

find_package(Glog REQUIRED)
include_directories (${GLOG_INCLUDE_DIRS})
add_executable(aliaslda ${SOURCE_FILES})
target_link_libraries(aliaslda ${JSONCPP_LIBRARIES} ${Boost_LIBRARIES} ${GLOG_LIBRARIES})