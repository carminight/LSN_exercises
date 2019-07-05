CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o DataBlocking.o funcvirtual.o Metropolis.o
	$(CC) random.o main.o DataBlocking.o funcvirtual.o Metropolis.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
DataBlocking.o : DataBlocking.cpp DataBlocking.h
	$(CC) -c DataBlocking.cpp -o DataBlocking.o $(CFLAGS)
funcvirtual.o : funcvirtual.cpp funcvirtual.h
	$(CC) -c funcvirtual.cpp -o funcvirtual.o $(CFLAGS)
Metropolis.o : Metropolis.cpp Metropolis.h
	$(CC) -c Metropolis.cpp -o Metropolis.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out *.out

run:
	./clean.sh && make && ./main.exe 
