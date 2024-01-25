#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.21)
project(SNneutrinos)
#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
#----------------------------------------------------------------------------
# Find ROOT variables if the variable G4ANALYSIS_USE_ROOT is set
#
find_package(ROOT REQUIRED)

if(${ROOT_FOUND})
	message ($ENV{ROOTSYS})
	message(STATUS "ROOT was found!")
endif(${ROOT_FOUND})




#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
message(${Geant4_INCLUDE_DIR})
message(${ROOT_INCLUDE_DIRS})
include_directories(${ROOT_INCLUDE_DIRS}
                    ${PROJECT_SOURCE_DIR}/include
		    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
file(GLOB macros ${PROJECT_SOURCE_DIR}/*.mac)
file(GLOB inputs ${PROJECT_SOURCE_DIR}/*.in)
#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(SNneutrinos SNneutrinos.cc ${sources} ${headers})
target_link_libraries(SNneutrinos ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build OpNovice. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(SNneutrinos_SCRIPTS
    macro.mac
    gui.mac
    vis.mac
  )

foreach(_script ${SNneutrinos_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS SNneutrinos DESTINATION bin)
install(FILES ${macros} ${inputs} DESTINATION bin)
