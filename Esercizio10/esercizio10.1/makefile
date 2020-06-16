CC = g++
CFLAGS = -Wall -O3 --std=c++11

es10.1.x : es10.1.o random.o GA.o
	$(CC) random.o GA.o es10.1.o -o es10.1.x
es10.1.o : es10.1.cc random.h GA.h
	$(CC) -c es10.1.cc -o es10.1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
GA.o : GA.cc GA.h random.h
	$(CC) -c GA.cc -o GA.o $(CFLAGS)

clean:
	rm *.o main.exe seed.out
