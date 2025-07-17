# src/alglib/CMakeLists.txt
add_library(alglib STATIC
    source/alglibinternal.cpp
    source/alglibmisc.cpp
    source/ap.cpp
    source/integration.cpp
    source/interpolation.cpp
    source/linalg.cpp
    source/optimization.cpp
    source/solvers.cpp
    source/specialfunctions.cpp
)
target_compile_options(alglib PRIVATE -w) # Disabilita warning per alglib

