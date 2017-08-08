CC = g++ -g
CC = g++ 

CFLAGS = -g -fpermissive -I..

.SUFFIXES: .cpp

.cpp.o:
	$(CC) $(CFLAGS) -c $*.cpp -o $*.o

OBJS = argprocessorTests.o autocorrTests.o basicTests.o clkTests.o randomTests.o unitTestMain.o recomboCriteriaTests.o

unitTests: $(OBJS)
	@echo Build unit tests.
	$(CC) -o ../bin/test $(CFLAGS) $(OBJS) ../librecombo.a ../libunittest.a

samples: sampleUnitTestMain.o sampleUnitTests.o
	@echo Build example tests.
	$(CC) -o ../bin/samples $(CFLAGS) sampleUnitTestMain.o sampleUnitTests.o ../libunittest.a

