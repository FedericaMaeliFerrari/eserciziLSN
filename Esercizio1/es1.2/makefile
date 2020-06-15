CC = g++
CFLAGS = -Wall -O3 --std=c++11

es1.2.x : es1.2.o random.o randomDistrib.o
	$(CC) random.o randomDistrib.o es1.2.o -o es1.2.x
es1.2.o : es1.2.cc random.h randomDistrib.h
	$(CC) -c es1.2.cc -o es1.2.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
randomDistrib.o : randomDistrib.cc randomDistrib.h
	$(CC) -c randomDistrib.cc -o randomDistrib.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
