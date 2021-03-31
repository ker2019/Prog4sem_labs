cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_FLAGS "-std=c++11 -g")

project(fenicsdemo)

find_package(DOLFIN REQUIRED)
include(${DOLFIN_USE_FILE})

#add_library(gmsh SHARED IMPORTED)
#set_target_properties(gmsh PROPERTIES
#   IMPORTED_LOCATION "/home/sandbox/.local/lib/python3.8/site-packages/gmsh-git-Linux64-sdk/lib/libgmsh.so"
#   INTERFACE_INCLUDE_DIRECTORIES "/home/sandbox/.local/lib/python3.8/site-packages/gmsh-git-Linux64-sdk/include"
#)
#set(GMSH_LIBRARY gmsh)

add_executable(wave wave.cpp)

target_link_libraries(wave dolfin)
#target_link_libraries(Difraction ${GMSH_LIBRARY})
