cmake_minimum_required(VERSION 3.25)

project(IntaRNAnew VERSION 4.0.0 LANGUAGES CXX)

set(CMAKE_CXX_EXTENSIONS OFF)

include(CTest)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
find_package(Threads REQUIRED)

option(INTARNANEW_ENABLE_SANITIZERS
    "Enable address and undefined-behaviour sanitizers" OFF)
option(INTARNANEW_BUILD_BENCHMARKS
    "Build the POSIX benchmark when BUILD_TESTING is disabled" OFF)
set(INTARNA_LEGACY_BIN "$ENV{INTARNA_LEGACY_BIN}" CACHE FILEPATH
    "Optional IntaRNA Legacy executable for the 17-case black-box parity test")

set(_intarnanew_local_parameter_dir
    "${CMAKE_CURRENT_SOURCE_DIR}/../.conda-env/share/ViennaRNA")
if(NOT EXISTS "${_intarnanew_local_parameter_dir}")
    set(_intarnanew_local_parameter_dir "")
endif()
set(INTARNANEW_TEST_PARAMETER_DIR "${_intarnanew_local_parameter_dir}" CACHE PATH
    "Directory containing public ViennaRNA .par files for thermodynamic tests")
unset(_intarnanew_local_parameter_dir)

function(intarnanew_enable_warnings target)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4 /permissive-)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
        target_compile_options(${target} PRIVATE
            -Wall -Wextra -Wpedantic -Wconversion -Wshadow)
    endif()
endfunction()

function(intarnanew_enable_sanitizers target)
    if(NOT INTARNANEW_ENABLE_SANITIZERS)
        return()
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
        target_compile_options(${target} PRIVATE
            -fsanitize=address,undefined -fno-omit-frame-pointer)
        target_link_options(${target} PRIVATE -fsanitize=address,undefined)
    else()
        message(FATAL_ERROR
            "INTARNANEW_ENABLE_SANITIZERS requires GCC, Clang, or AppleClang")
    endif()
endfunction()

function(intarnanew_set_parameter_test_environment test_name)
    if(INTARNANEW_TEST_PARAMETER_DIR)
        set_tests_properties(${test_name} PROPERTIES ENVIRONMENT
            "INTARNANEW_PARAMETER_DIR=${INTARNANEW_TEST_PARAMETER_DIR}")
    endif()
endfunction()

add_library(intarnanew STATIC
    src/accessibility.cpp
    src/cli.cpp
    src/compression.cpp
    src/energy.cpp
    src/folding.cpp
    src/folding_parameters.cpp
    src/helix_blocks.cpp
    src/output.cpp
    src/output_plan.cpp
    src/predictor.cpp
    src/runner.cpp
    src/sequence.cpp
    src/thermo_parameters.cpp
)
add_library(IntaRNAnew::core ALIAS intarnanew)
set_target_properties(intarnanew PROPERTIES EXPORT_NAME core)
target_include_directories(intarnanew PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_features(intarnanew PUBLIC cxx_std_23)
target_link_libraries(intarnanew PUBLIC Threads::Threads)
intarnanew_enable_warnings(intarnanew)
intarnanew_enable_sanitizers(intarnanew)

add_library(intarnanew_tools STATIC
    src/tools/csv.cpp
    src/tools/mutations.cpp
    src/tools/pvalue.cpp
    src/tools/statistics.cpp
    src/tools/svg.cpp
)
add_library(IntaRNAnew::tools ALIAS intarnanew_tools)
set_target_properties(intarnanew_tools PROPERTIES
    EXPORT_NAME tools
    OUTPUT_NAME intarnanew-tools
)
target_include_directories(intarnanew_tools PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_features(intarnanew_tools PUBLIC cxx_std_23)
target_link_libraries(intarnanew_tools PUBLIC intarnanew Threads::Threads)
intarnanew_enable_warnings(intarnanew_tools)
intarnanew_enable_sanitizers(intarnanew_tools)

add_executable(IntaRNAnew app/main.cpp)
target_link_libraries(IntaRNAnew PRIVATE intarnanew)
intarnanew_enable_warnings(IntaRNAnew)
intarnanew_enable_sanitizers(IntaRNAnew)

foreach(tool IN ITEMS csv stats svg mutate pvalue)
    add_executable(intarnanew-${tool} tools/intarnanew-${tool}.cpp)
    target_link_libraries(intarnanew-${tool} PRIVATE intarnanew_tools)
    intarnanew_enable_warnings(intarnanew-${tool})
    intarnanew_enable_sanitizers(intarnanew-${tool})
endforeach()

if(INTARNANEW_BUILD_BENCHMARKS AND NOT UNIX)
    message(FATAL_ERROR "The black-box benchmark requires POSIX process APIs")
endif()
if(UNIX AND (INTARNANEW_BUILD_BENCHMARKS OR BUILD_TESTING))
    add_executable(IntaRNAnewBenchmark benchmarks/compare.cpp)
    target_compile_features(IntaRNAnewBenchmark PRIVATE cxx_std_23)
    intarnanew_enable_warnings(IntaRNAnewBenchmark)
    intarnanew_enable_sanitizers(IntaRNAnewBenchmark)
endif()

if(BUILD_TESTING)
    add_executable(IntaRNAnewIntegrationTests tests/integration.cpp)
    target_link_libraries(IntaRNAnewIntegrationTests PRIVATE intarnanew)
    intarnanew_enable_warnings(IntaRNAnewIntegrationTests)
    intarnanew_enable_sanitizers(IntaRNAnewIntegrationTests)
    add_test(NAME real_interactions COMMAND IntaRNAnewIntegrationTests)
    intarnanew_set_parameter_test_environment(real_interactions)

    add_executable(IntaRNAnewCliIoOutputTests tests/cli_io_output.cpp)
    target_link_libraries(IntaRNAnewCliIoOutputTests PRIVATE intarnanew)
    intarnanew_enable_warnings(IntaRNAnewCliIoOutputTests)
    intarnanew_enable_sanitizers(IntaRNAnewCliIoOutputTests)
    add_test(NAME cli_io_output_regressions COMMAND IntaRNAnewCliIoOutputTests)

    add_executable(IntaRNAnewConfigRegistryTests tests/config_registry.cpp)
    target_link_libraries(IntaRNAnewConfigRegistryTests PRIVATE intarnanew)
    intarnanew_enable_warnings(IntaRNAnewConfigRegistryTests)
    intarnanew_enable_sanitizers(IntaRNAnewConfigRegistryTests)
    add_test(NAME config_registry COMMAND IntaRNAnewConfigRegistryTests)
    set_tests_properties(config_registry PROPERTIES ENVIRONMENT
        "INTARNANEW_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}")

    add_executable(IntaRNAnewAccessibilityOracleTests tests/accessibility_oracle.cpp)
    target_link_libraries(IntaRNAnewAccessibilityOracleTests PRIVATE intarnanew)
    intarnanew_enable_warnings(IntaRNAnewAccessibilityOracleTests)
    intarnanew_enable_sanitizers(IntaRNAnewAccessibilityOracleTests)
    add_test(NAME accessibility_oracle COMMAND IntaRNAnewAccessibilityOracleTests)

    add_executable(IntaRNAnewFoldingOracleTests tests/folding_oracle.cpp)
    target_link_libraries(IntaRNAnewFoldingOracleTests PRIVATE intarnanew)
    intarnanew_enable_warnings(IntaRNAnewFoldingOracleTests)
    intarnanew_enable_sanitizers(IntaRNAnewFoldingOracleTests)
    add_test(NAME folding_oracle COMMAND IntaRNAnewFoldingOracleTests)
    intarnanew_set_parameter_test_environment(folding_oracle)

    add_executable(IntaRNAnewPredictorOracleTests tests/predictor_oracle.cpp)
    target_link_libraries(IntaRNAnewPredictorOracleTests PRIVATE intarnanew)
    intarnanew_enable_warnings(IntaRNAnewPredictorOracleTests)
    intarnanew_enable_sanitizers(IntaRNAnewPredictorOracleTests)
    add_test(NAME predictor_oracle COMMAND IntaRNAnewPredictorOracleTests)

    add_executable(IntaRNAnewExecutionRegionsTests tests/execution_regions.cpp)
    target_link_libraries(IntaRNAnewExecutionRegionsTests PRIVATE intarnanew)
    intarnanew_enable_warnings(IntaRNAnewExecutionRegionsTests)
    intarnanew_enable_sanitizers(IntaRNAnewExecutionRegionsTests)
    add_test(NAME region_window_execution COMMAND IntaRNAnewExecutionRegionsTests)

    add_executable(IntaRNAnewModelSemanticsTests tests/model_semantics.cpp)
    target_link_libraries(IntaRNAnewModelSemanticsTests PRIVATE intarnanew)
    intarnanew_enable_warnings(IntaRNAnewModelSemanticsTests)
    intarnanew_enable_sanitizers(IntaRNAnewModelSemanticsTests)
    add_test(NAME interaction_model_semantics COMMAND IntaRNAnewModelSemanticsTests)

    add_executable(IntaRNAnewCompressionIoTests tests/compression_io.cpp)
    target_link_libraries(IntaRNAnewCompressionIoTests PRIVATE intarnanew)
    intarnanew_enable_warnings(IntaRNAnewCompressionIoTests)
    intarnanew_enable_sanitizers(IntaRNAnewCompressionIoTests)
    add_test(NAME compression_and_gzip_io COMMAND IntaRNAnewCompressionIoTests)

    add_executable(IntaRNAnewOutputExecutionTests tests/output_execution.cpp)
    target_link_libraries(IntaRNAnewOutputExecutionTests PRIVATE intarnanew)
    intarnanew_enable_warnings(IntaRNAnewOutputExecutionTests)
    intarnanew_enable_sanitizers(IntaRNAnewOutputExecutionTests)
    add_test(NAME output_planning_and_parallel_execution
        COMMAND IntaRNAnewOutputExecutionTests)

    add_executable(IntaRNAnewRunnerContractTests tests/runner_contracts.cpp)
    target_link_libraries(IntaRNAnewRunnerContractTests PRIVATE intarnanew)
    intarnanew_enable_warnings(IntaRNAnewRunnerContractTests)
    intarnanew_enable_sanitizers(IntaRNAnewRunnerContractTests)
    add_test(NAME in_process_pair_runner COMMAND IntaRNAnewRunnerContractTests)

    add_executable(IntaRNAnewToolsContractTests tests/tools_contracts.cpp)
    target_link_libraries(IntaRNAnewToolsContractTests PRIVATE intarnanew_tools)
    intarnanew_enable_warnings(IntaRNAnewToolsContractTests)
    intarnanew_enable_sanitizers(IntaRNAnewToolsContractTests)
    add_test(NAME native_tools_contracts COMMAND IntaRNAnewToolsContractTests)

    if(UNIX)
        add_executable(IntaRNAnewPvalueExecutableTests tests/pvalue_executable.cpp)
        target_compile_features(IntaRNAnewPvalueExecutableTests PRIVATE cxx_std_23)
        intarnanew_enable_warnings(IntaRNAnewPvalueExecutableTests)
        intarnanew_enable_sanitizers(IntaRNAnewPvalueExecutableTests)
        add_test(NAME pvalue_executable_contracts
            COMMAND IntaRNAnewPvalueExecutableTests
                "$<TARGET_FILE:intarnanew-pvalue>")
        set_tests_properties(pvalue_executable_contracts PROPERTIES TIMEOUT 120)

        set(_intarnanew_legacy_test_arguments)
        if(INTARNA_LEGACY_BIN)
            list(APPEND _intarnanew_legacy_test_arguments
                --legacy "${INTARNA_LEGACY_BIN}")
        endif()
        add_test(NAME legacy_blackbox_17
            COMMAND IntaRNAnewBenchmark
                --new "$<TARGET_FILE:IntaRNAnew>"
                --cases "${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/parameters"
                --verify-only --expect-cases 17 --skip-without-legacy
                ${_intarnanew_legacy_test_arguments})
        set_tests_properties(legacy_blackbox_17 PROPERTIES
            SKIP_REGULAR_EXPRESSION "^SKIP:")
        unset(_intarnanew_legacy_test_arguments)
    endif()
endif()

install(TARGETS intarnanew intarnanew_tools
    EXPORT IntaRNAnewTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
install(TARGETS IntaRNAnew
    intarnanew-csv intarnanew-stats intarnanew-svg intarnanew-mutate
    intarnanew-pvalue
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(DIRECTORY include/intarnanew
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.hpp"
)
install(FILES README.md MATHEMATICAL_AUDIT.md LICENSE
    DESTINATION ${CMAKE_INSTALL_DOCDIR}
)
install(FILES tools/README.md
    DESTINATION ${CMAKE_INSTALL_DOCDIR}/tools
)

set(intarnanew_package_dir
    "${CMAKE_INSTALL_LIBDIR}/cmake/IntaRNAnew")
configure_package_config_file(
    IntaRNAnewConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/IntaRNAnewConfig.cmake"
    INSTALL_DESTINATION "${intarnanew_package_dir}"
)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/IntaRNAnewConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)
install(EXPORT IntaRNAnewTargets
    FILE IntaRNAnewTargets.cmake
    NAMESPACE IntaRNAnew::
    DESTINATION "${intarnanew_package_dir}"
)
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/IntaRNAnewConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/IntaRNAnewConfigVersion.cmake"
    DESTINATION "${intarnanew_package_dir}"
)
