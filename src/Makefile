CC = g++ -g
CC = g++ 

CFLAGS = -g -fpermissive

.SUFFIXES: .cpp

.cpp.o:
	$(CC) $(CFLAGS) -c $*.cpp -o $*.o

OBJS = autocorr.o bfacfProbabilities.o bfacfProbabilitiesFromZ.o bfacfProbabilitiesFromZFixed.o clk.o clkCigar.o clkConformationAsList.o clkConformationBfacf3.o clkConformationSeqconvert.o conformationAsList.o conformationWlc.o conformation.o legacyBfacf.o legacy.o mmchain.o random.o recomboFromFile.o regularNgon.o 

TESTOBJS = testFramework.o unitTests.o unitTestMain.o

all: programs test 

programs: mmc recomboFromFile

mmc: $(OBJS) mmcMain.o
	@echo Build mmc.
	$(CC) -o mmc $(CFLAGS) $(OBJS) mmcMain.o

recomboFromFile: $(OBJS) recomboFromFileMain.o
	@echo Build recombo from file.
	$(CC) -o recomboFromFile  $(CFLAGS) $(OBJS) recomboFromFileMain.o

unitTests: $(OBJS) $(TESTOBJS)
	@echo Build unit tests.
	$(CC) -o unitTests $(CFLAGS) $(OBJS) $(TESTOBJS)

test: unitTests
	@echo Run unit tests.
	./unitTests

clean: 
	@echo Clean.
	rm -f *.o

clobber: clean
	@echo Clobber.
	rm -f mmc recomboFromFile unitTests

