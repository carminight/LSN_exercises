CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o DataBlocking.o Buffon.o
	$(CC) random.o main.o DataBlocking.o Buffon.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
DataBlocking.o : DataBlocking.cpp DataBlocking.h
	$(CC) -c DataBlocking.cpp -o DataBlocking.o $(CFLAGS)
Buffon.o : Buffon.cpp Buffon.h
	$(CC) -c Buffon.cpp -o Buffon.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
run:
	make && ./main.exe