# Compile everything under 'dependencies'
add_subdirectory(deps)

# Add include directories
include_directories(
	"." 
	deps
	deps/yaml-cpp/include
)

# Find source files
aux_source_directory(
	"." 
	waternet-src
)

aux_source_directory(
	"deps" 
	deps_src
)

# Add compilation target
add_executable( 
	waternet
	${waternet-src} 
	${deps_src}
)

add_library(waternet-lib
        SHARED
	${waternet-src} 

	${deps_src}
)

target_link_libraries (
	waternet 
	yaml-cpp
)

target_link_libraries (
	waternet-lib
	yaml-cpp
)