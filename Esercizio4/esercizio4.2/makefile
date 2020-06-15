CC = g++
CFLAGS = -Wall -O3 --std=c++11

es4.2.x : es4.2.o MD.o random.o 
	$(CC) random.o MD.o es4.2.o -o es4.2.x
es4.2.o : es4.2.cc MD.h random.h
	$(CC) -c es4.2.cc -o es4.2.o $(CFLAGS)
MD.o : MD.cc MD.h random.h 
	$(CC) -c MD.cc -o MD.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean:
	rm -rf *epot*
	rm -rf *ekin*
	rm -rf *etot*
	rm -rf *t*

