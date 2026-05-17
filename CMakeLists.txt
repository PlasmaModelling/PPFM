cmake_minimum_required(VERSION 3.16)
project(PPFM LANGUAGES CXX)

# === C++ standard ===
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# === Use ccache if available ===
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
endif()

# === Enable FetchContent ===
include(FetchContent)

# === Detect available cores ===
include(ProcessorCount)
ProcessorCount(N)
if(NOT N EQUAL 0)
    set(CMAKE_BUILD_PARALLEL_LEVEL ${N})
endif()

# === Python ===
find_package(Python3 REQUIRED COMPONENTS Interpreter)

# ============================================================
# Main selection
# Default: main.cpp
# User may override with:
#   cmake -B build -S . -DMAIN=demo/mainArH2.cpp
# ============================================================
# ============================================================
# Main source selection
# ============================================================
if(DEFINED CACHE{MAIN})
    set(MAIN "$CACHE{MAIN}")
else()
    set(MAIN "main.cpp")
endif()

unset(MAIN CACHE)

# Resolve absolute path
if(IS_ABSOLUTE "${MAIN}")
    set(SELECTED_MAIN_SOURCE "${MAIN}")
else()
    set(SELECTED_MAIN_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/${MAIN}")
endif()

if(NOT EXISTS "${SELECTED_MAIN_SOURCE}")
    message(FATAL_ERROR "Main source file not found: ${SELECTED_MAIN_SOURCE}")
endif()

get_filename_component(SELECTED_MAIN_NAME_WE "${SELECTED_MAIN_SOURCE}" NAME_WE)
set(SELECTED_TARGET_NAME "${SELECTED_MAIN_NAME_WE}")
set(SELECTED_OUTPUT_NAME "${SELECTED_MAIN_NAME_WE}.out")

message(STATUS "PPFM main source selected: ${SELECTED_MAIN_SOURCE}")
message(STATUS "Target name: ${SELECTED_TARGET_NAME}")
message(STATUS "Output executable: ${SELECTED_OUTPUT_NAME}")
if(IS_ABSOLUTE "${MAIN}")
    set(SELECTED_MAIN_SOURCE "${MAIN}")
else()
    set(SELECTED_MAIN_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/${MAIN}")
endif()

if(NOT EXISTS "${SELECTED_MAIN_SOURCE}")
    message(FATAL_ERROR "Selected main source does not exist: ${SELECTED_MAIN_SOURCE}")
endif()

get_filename_component(SELECTED_MAIN_NAME_WE "${SELECTED_MAIN_SOURCE}" NAME_WE)

set(SELECTED_TARGET_NAME "${SELECTED_MAIN_NAME_WE}")
set(SELECTED_OUTPUT_NAME "${SELECTED_MAIN_NAME_WE}.out")
set(PPFM_EXECUTABLES_DIR "${CMAKE_CURRENT_SOURCE_DIR}/executables")

file(MAKE_DIRECTORY ${PPFM_EXECUTABLES_DIR})

message(STATUS "PPFM main source selected: ${SELECTED_MAIN_SOURCE}")
message(STATUS "Target name: ${SELECTED_TARGET_NAME}")
message(STATUS "Output executable: ${SELECTED_OUTPUT_NAME}")
message(STATUS "Executables output dir: ${PPFM_EXECUTABLES_DIR}")

# ============================================================
# Dynamically generate AcceptedSpecies.h from selected main
# ============================================================
set(ACCEPTED_SPECIES_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/src/AcceptedSpecies.h")
set(ACCEPTED_SPECIES_SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/src/AcceptedSpeciesGuard.py")

add_custom_target(generate_species_header ALL
    COMMAND ${CMAKE_COMMAND} -E echo "Generating AcceptedSpecies.h from ${SELECTED_MAIN_SOURCE}"
    COMMAND ${Python3_EXECUTABLE} ${ACCEPTED_SPECIES_SCRIPT} ${SELECTED_MAIN_SOURCE}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src
    COMMENT "Running AcceptedSpeciesGuard.py"
)

# ============================================================
# External library: alglib
# ============================================================
add_subdirectory(src/alglib)

# ============================================================
# Project sources
# ============================================================
file(GLOB_RECURSE PPFM_SOURCES CONFIGURE_DEPENDS "src/*.cpp")

# ============================================================
# PPFM core library
# ============================================================
add_library(ppfm_core ${PPFM_SOURCES})

target_compile_definitions(ppfm_core PUBLIC
    PPFM_PROJECT_ROOT=\"${CMAKE_SOURCE_DIR}\"
)

target_include_directories(ppfm_core PUBLIC ${CMAKE_SOURCE_DIR}/src)
add_dependencies(ppfm_core generate_species_header)
target_link_libraries(ppfm_core PUBLIC alglib)

# === MSVC: full C++20 support ===
if(MSVC)
    target_compile_options(ppfm_core PRIVATE /std:c++20 /Zc:__cplusplus)
endif()

# === OpenMP ===
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(ppfm_core PUBLIC OpenMP::OpenMP_CXX)
endif()

# === libcurl ===
if(UNIX AND NOT APPLE)
    find_package(CURL QUIET)
    if(NOT CURL_FOUND)
        message(STATUS "libcurl not found, downloading through FetchContent...")
        FetchContent_Declare(
            curl
            URL https://github.com/curl/curl/releases/download/curl-8_8_0/curl-8.8.0.zip
        )
        FetchContent_MakeAvailable(curl)
        set(CURL_FOUND TRUE)
    endif()
endif()

if(CURL_FOUND)
    target_link_libraries(ppfm_core PUBLIC CURL::libcurl)
    target_compile_definitions(ppfm_core PUBLIC PPFM_USE_CURL)
else()
    message(WARNING "libcurl not found. Network-dependent features will be disabled.")
endif()

# ============================================================
# Main executable
# ============================================================
add_executable(${SELECTED_TARGET_NAME} ${SELECTED_MAIN_SOURCE})
add_custom_command(TARGET ${SELECTED_TARGET_NAME}
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E touch $<TARGET_FILE:${SELECTED_TARGET_NAME}>
)
add_custom_target(force_${SELECTED_TARGET_NAME} ALL
    DEPENDS ${SELECTED_TARGET_NAME}
)
add_dependencies(force_${SELECTED_TARGET_NAME} ${SELECTED_TARGET_NAME})
target_link_libraries(${SELECTED_TARGET_NAME} PRIVATE ppfm_core)

set_target_properties(${SELECTED_TARGET_NAME} PROPERTIES
    OUTPUT_NAME "${SELECTED_OUTPUT_NAME}"
    RUNTIME_OUTPUT_DIRECTORY "${PPFM_EXECUTABLES_DIR}"
)

# === For MSVC: place output binary in executables/ for all configurations ===
if(MSVC)
    set_target_properties(${SELECTED_TARGET_NAME} PROPERTIES
        OUTPUT_NAME "${SELECTED_OUTPUT_NAME}"
        RUNTIME_OUTPUT_DIRECTORY_DEBUG "${PPFM_EXECUTABLES_DIR}"
        RUNTIME_OUTPUT_DIRECTORY_RELEASE "${PPFM_EXECUTABLES_DIR}"
        RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO "${PPFM_EXECUTABLES_DIR}"
        RUNTIME_OUTPUT_DIRECTORY_MINSIZEREL "${PPFM_EXECUTABLES_DIR}"
    )
endif()

# === Linux linker flags ===
if(UNIX AND NOT APPLE)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fuse-ld=lld")
endif()

option(PPFM_BUILD_TESTS "Build PPFM tests" OFF)

if(PPFM_BUILD_TESTS)
    include(CTest)
    add_subdirectory(tests)
endif()