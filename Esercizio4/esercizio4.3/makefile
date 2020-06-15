CC = g++
CFLAGS = -Wall -O3 --std=c++11

es4.3.x : es4.3.o MD.o random.o 
	$(CC) MD.o random.o es4.3.o -o es4.3.x
es4.3.o : es4.3.cc MD.h random.h
	$(CC) -c es4.3.cc -o es4.3.o $(CFLAGS)
MD.o : MD.cc MD.h random.h
	$(CC) -c MD.cc -o MD.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean:
	rm -rf *epot*
	rm -rf *ekin*
	rm -rf *etot*
	rm -rf *t*

