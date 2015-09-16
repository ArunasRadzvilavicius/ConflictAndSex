CC=g++
CCOPT=-O3 -Wall

all : Unisex
	echo `date` @ `uname -n` > VERSION

Unisex : Unisex.o Moose.o
	$(CC) $(CCOPT) -o unisex $^ -lgsl -lgslcblas -lm

Unisex.o : Unisex.cpp Moose.h
	$(CC) $(CCOPT) -c Unisex.cpp

Moose.o : Moose.cpp Moose.h
	$(CC) $(CCOPT) -c Moose.cpp

clean :
	rm unisex *.o
