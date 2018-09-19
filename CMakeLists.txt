cmake_minimum_required(VERSION 3.8)
project(recombo)
include_directories(src)
set(CMAKE_CXX_STANDARD 11)
set(EXECUTABLE_OUTPUT_PATH src/bin)
set(TEST_SOURCE_FILES
        src/tests/argprocessorTests.cpp
        src/tests/argprocessorTests.h
        src/tests/autocorrTests.cpp
        src/tests/autocorrTests.h
        src/tests/basicTests.cpp
        src/tests/basicTests.h
        src/tests/clkTests.cpp
        src/tests/clkTests.h
        src/tests/randomTests.cpp
        src/tests/randomTests.h
        src/tests/recomboCriteriaTests.cpp
        src/tests/recomboCriteriaTests.h
        src/tests/sampleUnitTests.cpp
        src/tests/sampleUnitTests.h
        src/tests/testFramework.cpp
        src/tests/testFramework.h)
set(SHARED_SOURCE_FILES
        src/argprocessor.cpp
        src/argprocessor.h
        src/autocorr.cpp
        src/autocorr.h
        src/bfacfProbabilities.cpp
        src/bfacfProbabilities.h
        src/bfacfProbabilitiesFromZ.cpp
        src/bfacfProbabilitiesFromZ.h
        src/bfacfProbabilitiesFromZFixed.cpp
        src/bfacfProbabilitiesFromZFixed.h
        src/clk.cpp
        src/clk.h
        src/clkCigar.cpp
        src/clkCigar.h
        src/clkConformationAsList.cpp
        src/clkConformationAsList.h
        src/clkConformationBfacf3.cpp
        src/clkConformationBfacf3.h
        src/clkConformationSeqconvert.cpp
        src/clkConformationSeqconvert.h
        src/conformation.cpp
        src/conformation.h
        src/conformationAsList.cpp
        src/conformationAsList.h
        src/conformationWlc.cpp
        src/conformationWlc.h
        src/genericConformation.h
        src/id.cpp
        src/id.h
        src/legacy.cpp
        src/legacy.h
        src/legacyBfacf.cpp
        src/legacyBfacf.h
        src/main.cpp
        src/mmchain.cpp
        src/mmchain.h
        src/polynomesAZ.c
        src/polynomesAZ.h
        src/pseudorandom.cpp
        src/pseudorandom.h
        src/recomboArgprocessor.cpp
        src/recomboArgprocessor.h
        src/recomboFromFile.cpp
        src/recomboFromFile.h
        src/regularNgon.cpp
        src/regularNgon.h
        src/runRecombo.cpp
        src/runRecombo.h
        src/threevector.h
        src/zAnalyzer.cpp
        src/zAnalyzer.h
        config.h)
set(XINGER_SOURCE_FILES
        src/xinger/random.cpp
        src/xinger/random.h
        src/xinger/util.cpp
        src/xinger/util.h
        src/xinger/xinger.cpp
        src/xinger/xinger.h)

add_executable(recomboFromFile ${SHARED_SOURCE_FILES} src/recomboFromFileMain.cpp)
add_executable(xinger ${XINGER_SOURCE_FILES} src/xinger/main.cpp)
add_executable(zAnalyzer ${SHARED_SOURCE_FILES} src/zAnalyzerMain.cpp)
add_executable(unitTest ${SHARED_SOURCE_FILES} ${TEST_SOURCE_FILES} src/tests/unitTestMain.cpp)
add_executable(mmc ${SHARED_SOURCE_FILES} src/mmcMain.cpp)
add_executable(homfly ${SHARED_SOURCE_FILES} src/homfly1_21.c src/cross.c src/cross.h)
add_executable(idknot ${SHARED_SOURCE_FILES} src/idknot.cpp)

