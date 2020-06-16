CC = g++
CFLAGS = -Wall -O3 --std=c++11

es7.2.x : es7.2.o random.o NVT.o
	$(CC) random.o NVT.o es7.2.o -o es7.2.x
es7.2.o : es7.2.cc random.h NVT.h
	$(CC) -c es7.2.cc -o es7.2.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
NVT.o : NVT.cc NVT.h random.h
	$(CC) -c NVT.cc -o NVT.o $(CFLAGS)

clean:
	rm *.o main.exe seed.out
