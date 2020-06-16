CC = g++
CFLAGS = -Wall -O3 --std=c++11

es10.2.x : es10.2.o random.o GA.o
	$(CC) random.o GA.o es10.2.o -o es10.2.x
es10.2.o : es10.2.cc random.h GA.h
	$(CC) -c es10.2.cc -o es10.2.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
GA.o : GA.cc GA.h random.h
	$(CC) -c GA.cc -o GA.o $(CFLAGS)

clean:
	rm *.o main.exe seed.out
