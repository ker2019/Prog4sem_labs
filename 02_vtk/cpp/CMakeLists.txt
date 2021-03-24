cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_FLAGS "-std=c++11 -g")

project(vtkdemo)
add_executable(heart heart.cpp)

set(VTK_DIR /home/sandbox/vtkbase/lib/cmake/vtk-7.1)
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_library(gmsh SHARED IMPORTED)
set_target_properties(gmsh PROPERTIES
    IMPORTED_LOCATION "/home/sandbox/.local/lib/python3.8/site-packages/gmsh-git-Linux64-sdk/lib/libgmsh.so"
    INTERFACE_INCLUDE_DIRECTORIES "/home/sandbox/.local/lib/python3.8/site-packages/gmsh-git-Linux64-sdk/include"
)
set(GMSH_LIBRARY gmsh)

target_link_libraries(heart ${VTK_LIBRARIES})
target_link_libraries(heart ${GMSH_LIBRARY})
