CC = g++
CFLAGS = -Wall -O3 --std=c++11

es7.1.x : es7.1.o random.o NVT.o
	$(CC) random.o NVT.o es7.1.o -o es7.1.x
es7.1.o : es7.1.cc random.h NVT.h
	$(CC) -c es7.1.cc -o es7.1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
NVT.o : NVT.cc NVT.h random.h
	$(CC) -c NVT.cc -o NVT.o $(CFLAGS)

clean:
	rm *.o main.exe seed.out
