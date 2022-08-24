
include(FetchContent)
FetchContent_Declare(
    catamari
    GIT_REPOSITORY https://gitlab.com/hodge_star/catamari.git
    GIT_TAG 47e6369cb7139640ec74cec74020f2b3ebdf5983
    GIT_SHALLOW FALSE
)

FetchContent_Declare(
    mantis
    GIT_REPOSITORY https://gitlab.com/hodge_star/mantis.git
    GIT_TAG ceaed6ffb766b1e6a3cdcac726455d01525d8cd9
    GIT_SHALLOW FALSE
)

FetchContent_Declare(
    quotient
    GIT_REPOSITORY https://gitlab.com/hodge_star/quotient.git
    GIT_TAG cf8a9ae5853821e62320bb0be6c87fa6ad1dec36
    GIT_SHALLOW FALSE
)

FetchContent_MakeAvailable(catamari mantis quotient)
