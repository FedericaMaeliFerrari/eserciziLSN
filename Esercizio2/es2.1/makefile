CC = g++
CFLAGS = -Wall -O3 --std=c++11

es2.1.x : es2.1.o random.o funzione.o integrale.o funzioneBase.h
	$(CC) random.o funzione.o integrale.o es2.1.o -o es2.1.x
es2.1.o : es2.1.cc random.h funzione.h integrale.h funzioneBase.h 
	$(CC) -c es2.1.cc -o es2.1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
funzione.o : funzione.cc funzione.h funzioneBase.h
	$(CC) -c funzione.cc -o funzione.o $(CFLAGS)
integrale.o : integrale.cc random.h funzioneBase.h integrale.h
	$(CC) -c integrale.cc -o integrale.o $(CFLAGS)
clean:
	rm *.o main.exe seed.out
