CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o DataBlocking.o Option.o
	$(CC) random.o main.o DataBlocking.o Option.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
DataBlocking.o : DataBlocking.cpp DataBlocking.h
	$(CC) -c DataBlocking.cpp -o DataBlocking.o $(CFLAGS)
Option.o : Option.cpp Option.h
	$(CC) -c Option.cpp -o Option.o $(CFLAGS)

clean :
	rm *.o main.exe seed.out
run :
	make && ./main.exe