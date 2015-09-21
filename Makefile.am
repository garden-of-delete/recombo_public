AM_CXXFLAGS = -g -fpermissive -Isrc
AM_CFLAGS =
AM_LDFLAGS =

SRCS = src/argprocessor.cpp src/autocorr.cpp src/bfacfProbabilities.cpp src/bfacfProbabilitiesFromZ.cpp src/bfacfProbabilitiesFromZFixed.cpp src/clk.cpp src/clkCigar.cpp src/clkConformationAsList.cpp src/clkConformationBfacf3.cpp src/clkConformationSeqconvert.cpp src/conformationAsList.cpp src/conformationWlc.cpp src/conformation.cpp src/legacyBfacf.cpp src/legacy.cpp src/mmchain.cpp src/pseudorandom.cpp src/recomboArgprocessor.cpp src/recomboFromFile.cpp src/regularNgon.cpp src/runRecombo.cpp src/zAnalyzer.cpp src/id.cpp

all: mmc recomboFromFile zAnalyzer idknot homfly xinger testrecombo test

bin_PROGRAMS = mmc recomboFromFile zAnalyzer idknot homfly xinger testrecombo
idknot_SOURCES = src/id.cpp src/idknot.cpp
homfly_SOURCES = src/cross.c src/polynomesAZ.c src/homfly1_21.c
mmc_SOURCES = $(SRCS) src/mmcMain.cpp
recomboFromFile_SOURCES = $(SRCS) src/recomboFromFileMain.cpp
zAnalyzer_SOURCES = $(SRCS) src/zAnalyzerMain.cpp

testrecombo_SOURCES = $(SRCS) src/testFramework.cpp src/tests/argprocessorTests.cpp src/tests/autocorrTests.cpp src/tests/basicTests.cpp src/tests/clkTests.cpp src/tests/randomTests.cpp src/tests/unitTestMain.cpp


xinger_SOURCES = src/xinger/main.cpp src/xinger/random.cpp src/xinger/util.cpp src/xinger/xinger.cpp
xinger_CXXFLAGS = $(AM_CXXFLAGS) -Isrc/xinger

ACLOCAL_AMFLAGS = -I m4 --install
EXTRA_DIST = m4/NOTES

test: testrecombo
	./testrecombo

