add_executable(test_arh2_light
    arh2_light.cpp
)

target_include_directories(test_arh2_light PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(test_arh2_light PRIVATE
    ppfm_core
)

add_test(
    NAME arh2_light
    COMMAND test_arh2_light
)