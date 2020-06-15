CC = g++
CFLAGS = -Wall -O3 --std=c++11

es5.1.x : es5.1.o random.o randomDistrib.o funzione.o vettore.o metropolis.o funzioneBase.h
	$(CC) random.o randomDistrib.o funzione.o vettore.o metropolis.o es5.1.o -o es5.1.x
es5.1.o : es5.1.cc random.h randomDistrib.h funzione.h vettore.h metropolis.h funzioneBase.h
	$(CC) -c es5.1.cc -o es5.1.o $(CFLAGS)
metropolis.o : metropolis.cc random.h randomDistrib.h funzioneBase.h
	$(CC) -c metropolis.cc -o metropolis.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
randomDistrib.o : randomDistrib.cc randomDistrib.h
	$(CC) -c randomDistrib.cc -o randomDistrib.o $(CFLAGS)
funzione.o : funzione.cc funzione.h funzioneBase.h
	$(CC) -c funzione.cc -o funzione.o $(CFLAGS)
vettore.o: vettore.cc vettore.h 
	$(CC) -c vettore.cc -o vettore.o $(CFLAGS)

clean:
	rm *.o main.exe seed.out
