CC = g++
CFLAGS = -Wall -O3 --std=c++11

es1.1_c.x : es1.1_c.o random.o 
	$(CC) random.o es1.1_c.o -o es1.1_c.x
es1.1_c.o : es1.1_c.cc random.h
	$(CC) -c es1.1_c.cc -o es1.1_c.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
