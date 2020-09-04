CC = g++ 
CR = gcc

CFLAGS = -g -fpermissive
CRFLAGS = 

.SUFFIXES: .cpp

.cpp.o:
	$(CC) $(CFLAGS) -c $*.cpp -o $*.o

.c.o:
	$(CR) $(CRFLAGS) -c $*.c -o $*.o

OBJS = argprocessor.o autocorr.o bfacfProbabilities.o bfacfProbabilitiesFromZ.o bfacfProbabilitiesFromZFixed.o clk.o clkCigar.o clkConformationAsList.o clkConformationBfacf3.o clkConformationSeqconvert.o conformationAsList.o conformationWlc.o conformation.o legacyBfacf.o legacy.o mmchain.o pseudorandom.o regularNgon.o zAnalyzer.o id.o

all: programs 

programs: mmc idknot homfly xinger

bin:
	mkdir -p bin

librecombo: $(OBJS)
	ar r librecombo.a $(OBJS)

mmc: bin $(OBJS) mmcMain.o
	@echo Build mmc.
	$(CC) -o bin/mmc $(CFLAGS) $(OBJS) mmcMain.o

zAnalyzer: bin $(OBJS) zAnalyzerMain.o
	@echo Build zAnalyzer
	$(CC) -o bin/zAnalyzer $(CFLAGS) $(OBJS) zAnalyzerMain.o

idknot: bin $(OBJS) idknot.o
	@echo Build idknot
	$(CC) -o bin/idknot $(CRFLAGS) $(OBJS) idknot.o

homfly: bin cross.o polynomesAZ.o homfly1_21.o 
	@echo Build homfly
	$(CR) -o bin/homfly $(CRFLAGS) cross.o polynomesAZ.o homfly1_21.o

xinger: bin
	@echo Build unit tests.
	(cd xinger; $(MAKE))

clean: 
	@echo Clean.
	rm -f *.o *.a tests/*.o xinger/*.o

clobber: clean
	@echo Clobber.
	rm -fr bin

