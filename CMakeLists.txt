cmake_minimum_required(VERSION 3.16)
project(PPFM LANGUAGES CXX)

# === C++ standard ===
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# === Usa ccache se disponibile ===
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
endif()

# === Abilita FetchContent ===
include(FetchContent)

# === Rileva automaticamente i core ===
include(ProcessorCount)
ProcessorCount(N)
if(NOT N EQUAL 0)
    set(CMAKE_BUILD_PARALLEL_LEVEL ${N})
endif()

# === Python cross-platform ===
find_package(Python3 REQUIRED COMPONENTS Interpreter)

# === Auto-generazione AcceptedSpecies.h ===
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/src/AcceptedSpecies.h
    COMMAND ${CMAKE_COMMAND} -E echo "⚙️  Generating AcceptedSpecies.h..."
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/src/AcceptedSpeciesGuard.py
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src
    COMMENT "Running AcceptedSpeciesGuard.py to update AcceptedSpecies.h"
)
add_custom_target(generate_species_header DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/AcceptedSpecies.h)

# === Aggiungi libreria alglib ===
add_subdirectory(src/alglib)

# === Sorgenti progetto ===
file(GLOB_RECURSE PPFM_SOURCES "src/*.cpp")
set(MAIN_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp")

# === Output dell'eseguibile nella root ===
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR})

add_executable(rooster ${PPFM_SOURCES} ${MAIN_SOURCE})

# For MSVC: place output binary in source directory (not build/Debug)
if(MSVC)
    set_target_properties(rooster PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY_DEBUG   ${CMAKE_SOURCE_DIR}
        RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_SOURCE_DIR}
        RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO ${CMAKE_SOURCE_DIR}
        RUNTIME_OUTPUT_DIRECTORY_MINSIZEREL ${CMAKE_SOURCE_DIR}
    )
endif()

target_include_directories(rooster PRIVATE ${CMAKE_SOURCE_DIR}/src)
add_dependencies(rooster generate_species_header)
target_link_libraries(rooster PRIVATE alglib)

# === MSVC: Abilita pieno supporto a C++20 (std::numbers, ecc.) ===
if(MSVC)
    target_compile_options(rooster PRIVATE /std:c++20 /Zc:__cplusplus)
endif()

# === OpenMP cross-platform ===
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(rooster PRIVATE OpenMP::OpenMP_CXX)
endif()

# === libcurl ===
if(UNIX AND NOT APPLE)
    find_package(CURL QUIET)
    if(NOT CURL_FOUND)
        message(STATUS "libcurl not found, downloading through FetchContent...")
        include(FetchContent)
        FetchContent_Declare(
            curl
            URL https://github.com/curl/curl/releases/download/curl-8_8_0/curl-8.8.0.zip
        )
        FetchContent_MakeAvailable(curl)
        set(CURL_FOUND TRUE)
    endif()
endif()

if(CURL_FOUND)
    target_link_libraries(rooster PRIVATE CURL::libcurl)
    target_compile_definitions(rooster PRIVATE PPFM_USE_CURL)
else()
    message(WARNING "!! libcurl not found. Network-dependent features will be disabled.")
endif()

# === Linker flags (solo Linux) ===
if(UNIX AND NOT APPLE)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fuse-ld=lld")
endif()
