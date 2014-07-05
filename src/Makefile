CC = g++ -g
CC = g++ 

CFLAGS = -g -fpermissive

.SUFFIXES: .cpp

.cpp.o:
	$(CC) $(CFLAGS) -c $*.cpp -o $*.o

OBJS = argprocessor.o autocorr.o bfacfProbabilities.o bfacfProbabilitiesFromZ.o bfacfProbabilitiesFromZFixed.o clk.o clkCigar.o clkConformationAsList.o clkConformationBfacf3.o clkConformationSeqconvert.o conformationAsList.o conformationWlc.o conformation.o legacyBfacf.o legacy.o mmchain.o pseudorandom.o recomboFromFile.o regularNgon.o zAnalyzer.o

all: programs librecombo libunittest test 

programs: mmc recomboFromFile zAnalyzer

bin:
	mkdir -p bin

librecombo: $(OBJS)
	ar r librecombo.a $(OBJS)

mmc: bin $(OBJS) mmcMain.o
	@echo Build mmc.
	$(CC) -o bin/mmc $(CFLAGS) $(OBJS) mmcMain.o

recomboFromFile: bin $(OBJS) recomboFromFileMain.o
	@echo Build recombo from file.
	$(CC) -o bin/recomboFromFile  $(CFLAGS) $(OBJS) recomboFromFileMain.o

zAnalyzer: bin $(OBJS) zAnalyzerMain.o
	@echo Build zAnalyzer
	$(CC) -o bin/zAnalyzer $(CFLAGS) $(OBJS) zAnalyzerMain.o

libunittest: bin testFramework.o
	ar r libunittest.a testFramework.o

unitTests: bin librecombo libunittest
	@echo Build unit tests.
	(cd tests; $(MAKE))

test: unitTests
	@echo Run unit tests.
	bin/test

sampleTest: unitTests
	@echo Run example unit tests.
	bin/samples

clean: 
	@echo Clean.
	rm -f *.o *.a tests/*.o

clobber: clean
	@echo Clobber.
	rm -fr bin

