cmake_minimum_required( VERSION 2.6 )

project(waternet)

aux_source_directory(
	"src" 
	waternet-src
)

add_executable (
	waternet 
	${waternet-src}
)

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
message("OPENMP FOUND")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_definitions( --std=c++11 )


target_link_libraries (
	waternet
	yaml-cpp
)