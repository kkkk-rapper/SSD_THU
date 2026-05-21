#----------------------------------------------------------------------------
# sim_project : HPGe gamma-spectroscopy simulation
#----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.16...3.27)
project(sim_project)

#----------------------------------------------------------------------------
# Find Geant4, activating all UI and Vis drivers by default.
#
find_package(Geant4 REQUIRED ui_all vis_all)

#----------------------------------------------------------------------------
# Sources and headers (headers listed so they show up in IDEs).
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(sim_project sim_project.cc ${sources} ${headers})
target_include_directories(sim_project PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(sim_project PRIVATE ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy macros next to the executable so it can be run from the build dir.
#
set(SIM_SCRIPTS
  init_vis.mac
  vis.mac
  run.mac
  )
foreach(_script ${SIM_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()
