CC = mpic++
CFLAGS = -Wall -O3 --std=c++17

main.exe : main.o random.o city.o
	$(CC) $(CFLAGS) city.o random.o main.o -o main.exe

main.o : main.cpp city.h random.h
	$(CC) $(CFLAGS) -c main.cpp -o main.o

random.o : random.cpp random.h
	$(CC) $(CFLAGS) -c random.cpp -o random.o

city.o : city.cpp city.h
	$(CC) $(CFLAGS) -c city.cpp -o city.o

clean :
	rm -f *.o main.exe seed.out
