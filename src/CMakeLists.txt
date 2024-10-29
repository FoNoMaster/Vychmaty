file(GLOB_RECURSE SOURCES *.hpp)
add_library(${PROJECT_NAME} INTERFACE ${SOURCES})

target_link_libraries(${PROJECT_NAME} INTERFACE Eigen3::Eigen)

target_include_directories(${PROJECT_NAME} INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_20)
