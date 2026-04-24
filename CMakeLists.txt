cmake_minimum_required(VERSION 3.16)
project(recombo)
set(CMAKE_CXX_STANDARD 11)
set(EXECUTABLE_OUTPUT_PATH src/bin)

include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.14.0
)
FetchContent_MakeAvailable(googletest)

# Core library shared by all executables and tests
add_library(recombo_core STATIC
        src/bfacf_config.cpp
        src/autocorr.cpp
        src/bfacfProbabilities.cpp
        src/bfacfProbabilitiesFromZ.cpp
        src/bfacfProbabilitiesFromZFixed.cpp
        src/clk.cpp
        src/clkCigar.cpp
        src/clkConformationAsList.cpp
        src/clkConformationBfacf3.cpp
        src/clkConformationSeqconvert.cpp
        src/conformation.cpp
        src/conformationAsList.cpp
        src/conformationWlc.cpp
        src/id.cpp
        src/json11.cpp
        src/legacy.cpp
        src/legacyBfacf.cpp
        src/mmc_config.cpp
        src/mmchain.cpp
        src/newsud.cpp
        src/polynomesAZ.c
        src/writhe.cpp
        src/pseudorandom.cpp
        src/recombo_config.cpp
        src/regularNgon.cpp
        src/zAnalyzer.cpp
        src/zanalyzer_config.cpp)
target_include_directories(recombo_core PUBLIC src)

# Executables
add_executable(xinger
        src/xinger/main.cpp
        src/xinger/random.cpp
        src/xinger/util.cpp
        src/xinger/xinger.cpp)
target_include_directories(xinger PRIVATE src/xinger)

add_executable(zAnalyzer src/zAnalyzerMain.cpp)
target_link_libraries(zAnalyzer recombo_core)

add_executable(mmc src/mmcMain.cpp)
target_link_libraries(mmc recombo_core)

add_executable(bfacf src/bfacfMain.cpp)
target_link_libraries(bfacf recombo_core)

add_executable(recombo src/recomboMain.cpp)
target_link_libraries(recombo recombo_core)

add_executable(homfly src/homfly1_21.c src/cross.c)
target_link_libraries(homfly recombo_core)
# Legacy K&R C files need relaxed compiler settings
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    set_source_files_properties(src/cross.c src/homfly1_21.c PROPERTIES
            COMPILE_OPTIONS "-Wno-implicit-int;-fcommon")
elseif(CMAKE_C_COMPILER_ID MATCHES "Clang")
    set_source_files_properties(src/cross.c src/homfly1_21.c PROPERTIES
            COMPILE_OPTIONS "-Wno-deprecated-non-prototype;-Wno-implicit-int")
endif()

add_executable(idknot src/idknot.cpp)
target_link_libraries(idknot recombo_core)

# Legacy unit tests
add_executable(legacyUnitTest
        src/tests/unitTestMain.cpp
        src/tests/autocorrTests.cpp
        src/tests/clkTests.cpp
        src/tests/randomTests.cpp
        src/tests/recomboCriteriaTests.cpp
        src/tests/testFramework.cpp)
target_link_libraries(legacyUnitTest recombo_core gtest gtest_main)
include(GoogleTest)
gtest_discover_tests(legacyUnitTest)

enable_testing()
add_subdirectory(test)
