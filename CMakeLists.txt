cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

# build on windows
# mkdir build
# cd build
# cmake -G "MinGW Makefiles" ..
# mingw32-make
# MaxwellConstruction.exe

# build on linux
# mkdir build
# cd build
# cmake ..
# make
# ./MaxwellConstruction

#=================== Set the example path ===========================
# Check the example directory to find more examples. By changing the 
# address below to the desired example, you set the compiler for it.
 
set(example_path ../examples/main.f90)
#====================================================================

#=================== Set project name ===============================
# The project executable file is made by this name.

set (project_name MaxwellConstruction)
#====================================================================


# The project name is applied.
project(${project_name})

# SET(CMAKE_C_COMPILER /path/to/c/compiler)
# SET(CMAKE_CXX_COMPILER /path/to/cpp/compiler)
# SET(CMAKE_Fortran_COMPILER /path/to/Fortran/compiler)

# Fortran language is activated.
enable_language(Fortran)

# You can set the build type in terminal as: 
#$ cmake -DCMAKE_BUILD_TYPE=Debug ..
# Otherwise "Release" is the default build type.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

# Gnu compiler flags
# -w                             switchs off all the warnings.
# -wall                          switchs on  all the warnings.
# -march=native                  optimization for the local machine
# -O0, -O1, -O2, -O3, -Ofast     overall optimization levels
# -g0,g1,g2,g3					 levels of debug information

# This flags are applied to all bluid types
set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -w")
# Debug flags
set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g3")
# Release flags
set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -march=native")

# Fortran module files are placed here.
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)


# All the subdirectories having CMakeLists.txt are listed.
add_subdirectory(src)


# The executable is built with the name of the project from the
# example path provided above.
add_executable(
    ${project_name}
    ${example_path}
    )

# Libraries needed by the executable are listed. 

target_link_libraries(
	${project_name} 
	src
)

# Display general status
message("---------------Status--------------")
message("Project name = ${project_name}")
message("Example path = ${example_path}")
message("Compiler = ${CMAKE_Fortran_COMPILER_ID}")
message("Build type = ${CMAKE_BUILD_TYPE}")
message("-----------------------------------")


