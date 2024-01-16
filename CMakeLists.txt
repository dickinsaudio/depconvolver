# CMakeList.txt : CMake project for DepConvolver, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

#SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -funsafe-math-optimizations -mfpu=neon")
#SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -funsafe-math-optimizations -D_GLIBCXX_USE_CXX11_ABI=0")

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -funsafe-math-optimizations -O3 -flto")

project ("DepConvolver")
add_subdirectory(daes67)

include_directories("." "~/intel/oneapi/ipp/latest/include")
# Add source to this project's executable.
add_executable (DepConvolver "DepConvolver.cpp" "ThreadedDSP.cpp")
target_link_libraries(DepConvolver daes67 "rt" "pthread")

#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

if (CMAKE_CXX_LIBRARY_ARCHITECTURE STREQUAL x86_64-linux-gnu)
    include_directories("/opt/intel/oneapi/ipp/latest/include/ipp")
    target_link_libraries(DepConvolver "/opt/intel/oneapi/ipp/latest/lib/libipps.a")
    target_link_libraries(DepConvolver "/opt/intel/oneapi/ipp/latest/lib/libippcore.a")
    target_link_libraries(DepConvolver "/opt/intel/oneapi/ipp/latest/lib/libippvm.a")
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "arm|aarch64")
    SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -funsafe-math-optimizations -mfpu=neon")
    target_link_libraries(DepConvolver "/home/ubuntu/git/Ne10/build/modules/libNE10.a")
endif()
set_property(TARGET DepConvolver PROPERTY CXX_STANDARD 17)
set_property(TARGET DepConvolver PROPERTY CXX_STANDARD_REQUIRED ON)
 
# TODO: Add tests and install targets if needed.
