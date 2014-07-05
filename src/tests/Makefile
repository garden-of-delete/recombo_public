CC = g++ -g
CC = g++ 

CFLAGS = -g -fpermissive -I..

.SUFFIXES: .cpp

.cpp.o:
	$(CC) $(CFLAGS) -c $*.cpp -o $*.o

OBJS = unitTestMain.o unitTests.o

all: samples 

samples: $(OBJS)
	@echo Build sample tests.
	$(CC) -o ../bin/samples $(CFLAGS) $(OBJS) ../libunittest.a

