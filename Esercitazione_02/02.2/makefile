CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o RW.o DataBlocking.o
	$(CC) random.o main.o RW.o DataBlocking.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
RW.o : RW.cpp RW.h
	$(CC) -c RW.cpp -o RW.o $(CFLAGS)
DataBlocking.o : DataBlocking.cpp DataBlocking.h
	$(CC) -c DataBlocking.cpp -o DataBlocking.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
run :
	make && ./main.exe