CC = g++ -g
CC = g++ 

# clean makefile for on-server solution build (rstolz 6/24/14)


CFLAGS = -g -fpermissive

.SUFFIXES: .cpp

.cpp.o:
	$(CC) $(CFLAGS) -c $*.cpp -o $*.o

OBJS = autocorr.o bfacfProbabilities.o bfacfProbabilitiesFromZ.o bfacfProbabilitiesFromZFixed.o clk.o clkCigar.o clkConformationAsList.o clkConformationBfacf3.o clkConformationSeqconvert.o conformationAsList.o conformationWlc.o conformation.o legacyBfacf.o legacy.o mmchain.o random.o recomboFromFile.o testFramework.o regularNgon.o 

mmc: 			$(OBJS) mmcMain.o
			echo build
			$(CC) -o mmc $(CFLAGS) $(OBJS) mmcMain.o

recomboFromFile: 	$(OBJS) recomboFromFileMain.o
			echo build
			$(cc) -o recomboFromFile  $(CFLAGS) $(OBJS) recomboFromFileMain.o

unitTests:		$(OBJS) unitTestMain.o
			echo build
			$(CC) -o unitTests $(CFLAGS) $(OBJS) unitTestMain.o

clean:
			rm $(OBJS) mmcMain.o unitTestMain.o recomboFromFileMain.o

