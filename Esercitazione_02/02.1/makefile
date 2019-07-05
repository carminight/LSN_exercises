CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o DataBlocking.o funcvirtual.o Integrate.o
	$(CC) random.o main.o DataBlocking.o funcvirtual.o  Integrate.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
DataBlocking.o : DataBlocking.cpp DataBlocking.h
	$(CC) -c DataBlocking.cpp -o DataBlocking.o $(CFLAGS)
funcvirtual.o : funcvirtual.cpp funcvirtual.h
	$(CC) -c funcvirtual.cpp -o funcvirtual.o $(CFLAGS)
Integrate.o : Integrate.cpp Integrate.h
	$(CC) -c Integrate.cpp -o Integrate.o $(CFLAGS)

clean :
	rm *.o main.exe seed.out
run:
	make && ./main.exe