CC = g++
CFLAGS = -Wall -O3 --std=c++11

es6.1.x : es6.1.o random.o Ising.o
	$(CC) random.o Ising.o es6.1.o -o es6.1.x
es6.1.o : es6.1.cc random.h Ising.h
	$(CC) -c es6.1.cc -o es6.1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
Ising.o : Ising.cc Ising.h random.h
	$(CC) -c Ising.cc -o Ising.o $(CFLAGS)

clean:
	rm *.o main.exe seed.out
