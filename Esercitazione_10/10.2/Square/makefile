CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o DataBlocking.o TravelingSalesman.o
	$(CC) random.o main.o DataBlocking.o TravelingSalesman.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
DataBlocking.o : DataBlocking.cpp DataBlocking.h
	$(CC) -c DataBlocking.cpp -o DataBlocking.o $(CFLAGS)
TravelingSalesman.o : TravelingSalesman.cpp TravelingSalesman.h
	$(CC) -c TravelingSalesman.cpp -o TravelingSalesman.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out *.out

run:
	./clean.sh && make && mpiexec -np 4 ./main.exe 
