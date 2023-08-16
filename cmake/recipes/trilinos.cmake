if(TARGET Trilinos::all_libs)
    return()
endif()

if(BUILD_TRILINOS_FROM_SOURCE)
    # Install Trilinos from source
    message(STATUS "Third-party: creating target 'Trilinos::all_libs'")

    set(Trilinos_ENABLE_Fortran OFF CACHE INTERNAL "" FORCE)
    set(TPL_ENABLE_MPI          OFF CACHE INTERNAL "" FORCE)
    set(Trilinos_ENABLE_OpenMP  ON CACHE INTERNAL "" FORCE)
    set(Trilinos_ENABLE_Belos   ON CACHE INTERNAL "" FORCE)
    set(Trilinos_ENABLE_Tpetra  ON CACHE INTERNAL "" FORCE)
    set(Trilinos_ENABLE_ML      ON CACHE INTERNAL "" FORCE)
    set(Trilinos_ENABLE_INSTALL_CMAKE_CONFIG_FILES ON CACHE INTERNAL "" FORCE)


    include(FetchContent)
    FetchContent_Declare(
        Trilinos
        GIT_REPOSITORY https://github.com/trilinos/Trilinos.git
        GIT_TAG trilinos-release-13-4-1
        GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(Trilinos)
    message(STATUS "Third-party: download Trilinos")
    add_library(Trilinos_Trilinos INTERFACE)
    add_library(Trilinos::Trilinos ALIAS Trilinos_Trilinos)
    target_include_directories(Trilinos_Trilinos INTERFACE ${Trilinos_INCLUDE_DIRS})
    target_link_libraries(Trilinos_Trilinos INTERFACE ${Trilinos_LIBRARIES})
else()
    find_package(Trilinos REQUIRED COMPONENTS ML Epetra)
    if(Trilinos_FOUND)
        message(STATUS "Third-party: found Trilinos")
        add_library(Trilinos_Trilinos INTERFACE)
        add_library(Trilinos::Trilinos ALIAS Trilinos_Trilinos)
        target_include_directories(Trilinos_Trilinos INTERFACE ${Trilinos_INCLUDE_DIRS})
        target_link_libraries(Trilinos_Trilinos INTERFACE ${Trilinos_LIBRARIES})
    else()
        message(STATUS "Third-party: Trilinos not found, disable Trilinos. Please build Trilinos from source or install it using your package manager.")
    endif()
endif()


