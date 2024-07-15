CC = g++
CFLAGS = -Wall -O3 --std=c++17

main.exe : main.o random.o city.o
	$(CC) city.o random.o main.o -o main.exe

main.o : main.cpp city.h random.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

city.o : city.cpp city.h
	$(CC) -c city.cpp -o city.o $(CFLAGS)

clean :
	rm *.o main.exe seed.out

