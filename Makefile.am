AM_CXXFLAGS = -g -fpermissive -std=c++11
AM_CFLAGS =
AM_LDFLAGS =

all: mmc zAnalyzer idknot homfly xinger testRecombo sampleTests test 

bin_PROGRAMS = mmc zAnalyzer idknot homfly xinger testRecombo sampleTests  

obj_src = src/argprocessor.cpp src/autocorr.cpp src/bfacfProbabilities.cpp src/bfacfProbabilitiesFromZ.cpp src/bfacfProbabilitiesFromZFixed.cpp src/clk.cpp src/clkCigar.cpp src/clkConformationAsList.cpp src/clkConformationBfacf3.cpp src/clkConformationSeqconvert.cpp src/conformationAsList.cpp src/conformationWlc.cpp src/conformation.cpp src/legacyBfacf.cpp src/legacy.cpp src/mmchain.cpp src/pseudorandom.cpp src/regularNgon.cpp src/zAnalyzer.cpp src/id.cpp

mmc_SOURCES = $(obj_src) src/mmcMain.cpp

zAnalyzer_SOURCES = $(obj_src) src/zAnalyzerMain.cpp

idknot_SOURCES = $(obj_src) src/idknot.cpp

homfly_SOURCES = $(obj_src) src/cross.c src/polynomesAZ.c src/homfly1_21.c

xinger_SOURCES = src/xinger/main.cpp src/xinger/random.cpp src/xinger/util.cpp src/xinger/xinger.cpp

testRecombo_SOURCES = $(obj_src) src/tests/testFramework.cpp src/tests/argprocessorTests.cpp src/tests/autocorrTests.cpp src/tests/clkTests.cpp src/tests/randomTests.cpp src/tests/recomboCriteriaTests.cpp src/tests/unitTestMain.cpp
testRecombo_CXXFLAGS = $(AM_CXXFLAGS) -Isrc

sampleTests_SOURCES = $(obj_src) src/tests/testFramework.cpp src/tests/sampleUnitTestMain.cpp src/tests/sampleUnitTests.cpp 
sampleTests_CXXFLAGS = $(AM_CXXFLAGS) -Isrc

ACLOCAL_AMFLAGS = -I m4 --install
EXTRA_DIST = m4/NOTES

test: all 
	./testRecombo
