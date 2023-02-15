# HYPRE GNU Lesser General Public License

if(TARGET HYPRE::HYPRE)
    return()
endif()

message(STATUS "Third-party: creating target 'HYPRE::HYPRE'")

set(HYPRE_WITH_MPI      OFF CACHE INTERNAL "" FORCE)
set(HYPRE_PRINT_ERRORS  ON CACHE INTERNAL "" FORCE)
set(HYPRE_BIGINT        ON CACHE INTERNAL "" FORCE)
set(HYPRE_USING_FEI     ON CACHE INTERNAL "" FORCE)
set(HYPRE_USING_OPENMP  ON CACHE INTERNAL "" FORCE)
set(HYPRE_SHARED       OFF CACHE INTERNAL "" FORCE)

include(FetchContent)
FetchContent_Declare(
    hypre
    GIT_REPOSITORY https://github.com/hypre-space/hypre.git
    GIT_TAG v2.25.0
    GIT_SHALLOW TRUE
)

FetchContent_MakeAvailable(hypre)

include_directories("${hypre_SOURCE_DIR}/src/FEI_mv")
add_subdirectory("${hypre_SOURCE_DIR}/src" ${hypre_BINARY_DIR})
file(REMOVE "${hypre_SOURCE_DIR}/src/utilities/version")