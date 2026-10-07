cmake_minimum_required(VERSION 3.25)

project(IntaRNAnewConsumerSmoke LANGUAGES CXX)

find_package(IntaRNAnew 1 CONFIG REQUIRED)

add_executable(intarnanew-consumer main.cpp)
target_link_libraries(intarnanew-consumer PRIVATE IntaRNAnew::tools)
target_compile_features(intarnanew-consumer PRIVATE cxx_std_23)
set_target_properties(intarnanew-consumer PROPERTIES CXX_EXTENSIONS OFF)
