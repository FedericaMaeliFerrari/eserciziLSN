CC = g++
CFLAGS = -Wall -O3 --std=c++11

es1.1.x : es1.1.o random.o 
	$(CC) random.o es1.1.o -o es1.1.x
es1.1.o : es1.1.cc random.h
	$(CC) -c es1.1.cc -o es1.1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
