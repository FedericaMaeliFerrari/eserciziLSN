CC = g++
CFLAGS = -Wall -O3 --std=c++11

es8.2.x : es8.2.o random.o Quantum.o
	$(CC) random.o Quantum.o es8.2.o -o es8.2.x
es8.1.o : es8.2.cc random.h Quantum.h
	$(CC) -c es8.2.cc -o es8.2.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
Quantum.o : Quantum.cc Quantum.h random.h
	$(CC) -c Quantum.cc -o Quantum.o $(CFLAGS)

clean:
	rm *.o main.exe seed.out
