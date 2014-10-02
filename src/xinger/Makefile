CC = g++ -g
CC = g++ 

CFLAGS = -g -fpermissive -I..

.SUFFIXES: .cpp

.cpp.o:
	$(CC) $(CFLAGS) -c $*.cpp -o $*.o

OBJS = main.o random.o util.o xinger.o

xinger: $(OBJS)
	$(CC) -o ../bin/xinger $(CFLAGS) $(OBJS)

