set(ENABLE_HELLOSPGRID OFF CACHE BOOL "Enable hello SPGrid project.")

if(ENABLE_HELLOSPGRID)
    add_subdirectory(demo_heat)
endif()