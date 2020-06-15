CC = g++
CFLAGS = -Wall -O3 --std=c++11

es3.x : es3.o random.o funzione.o integrale.o finanza.o randomDistrib.o funzioneBase.h
	$(CC) random.o funzione.o integrale.o finanza.o randomDistrib.o es3.o -o es3.x
es3.o : es3.cc random.h funzione.h integrale.h finanza.h funzioneBase.h randomDistrib.h
	$(CC) -c es3.cc -o es3.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
funzione.o : funzione.cc funzione.h funzioneBase.h
	$(CC) -c funzione.cc -o funzione.o $(CFLAGS)
integrale.o : integrale.cc random.h funzioneBase.h integrale.h
	$(CC) -c integrale.cc -o integrale.o $(CFLAGS)
finanza.o: finanza.cc finanza.h integrale.h randomDistrib.h
	$(CC) -c finanza.cc -o finanza.o $(CFLAGS)
randomDistrib.o: randomDistrib.cc randomDistrib.h
	$(CC) -c randomDistrib.cc -o randomDistrib.o $(CFLAGS)
clean:
	rm *.o main.exe seed.out
