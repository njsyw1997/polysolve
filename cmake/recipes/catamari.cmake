if(TARGET catamari::catamari)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    catamari
    GIT_REPOSITORY https://github.com/jpanetta/catamari.git
    GIT_TAG 7bfe01b8864e965f9e7a87150e848ceaa56e9cb7
    GIT_SHALLOW FALSE
)
FetchContent_MakeAvailable(catamari)


