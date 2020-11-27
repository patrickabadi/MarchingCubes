#cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

set(main MarchingCubesTest)
add_executable(${main} 
${main}.cpp 
Common.h 
FrameSet.cpp 
FrameSet.h 
MarchingCubes.cpp 
MarchingCubes.h
PclDefines.h
Visualize.cpp
Visualize.h )
target_link_libraries(${main} ${PCL_LIBRARIES})
set_target_properties(${main} PROPERTIES VS_DEBUGGER_COMMAND_ARGUMENTS "${PROJECT_SOURCE_DIR}/test_data/")
