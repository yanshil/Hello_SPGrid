add_executable(demo_heat
  main.cpp
  ../Heat_Transfer_Driver.cpp
  ../Heat_Transfer_Example.cpp
  )


find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

target_link_libraries(demo_heat
  PRIVATE
  nova::dynamics
  nova::geometry
  nova::spgrid
  nova::tools
  )

## Specify compile flags: only shows errors
target_compile_options(demo_heat PRIVATE -Wfatal-errors)

file( RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR} )

install(TARGETS demo_heat RUNTIME DESTINATION ${REL})